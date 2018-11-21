#include <H5Cpp.h>

#include "MedicalImages2HED.h"
#include "Common.h"
#include <caffe/caffe.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "itkOpenCVImageBridge.h"

#include "caffe/proto/caffe.pb.h"
#include "caffe/util/db.hpp"
#include "caffe/util/hdf5.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/rng.hpp"

using namespace caffe;
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

class Detector {

    public:
        Detector(const std::string & DataFolder,
                const std::string & DeployProtoFile,
                const std::string & TrainedModelFile, 
                const std::string & ListFile,
                int iDeviceId,
                int iOutputBlobIndex);
        
        bool Detect();
        bool DetectOne(std::string & strPatientId, int z);

        void SetColorOutput(bool bColor){
            m_bColor=bColor;
        }
        void SetThreshold(const float fThreshold){
            m_fThreshold = fThreshold;
        }

        typedef itk::Image<short, 2> ImageType;
        typedef itk::Image<short, 3> ImageType3D;

        virtual ~Detector() { }

    private:
        std::shared_ptr<caffe::Net<float> > m_p_net;
        nih::mi2hed * m_p_mi2hed;
        std::vector<std::string> m_vPatientIds;
        int m_channels;
        std::string m_strOutputName;
        bool m_bColor;
        float m_fThreshold;
};

Detector::Detector(const string & DataFolder,
                    const std::string & DeployProtoFile,
                    const std::string & TrainedModelFile, 
                    const std::string & ListFile,
                    int iDeviceId,
                    int iOutputBlobIndex){
    
    m_p_mi2hed = new nih::mi2hed;

    m_p_mi2hed->LoadPatientIds(ListFile);

    m_p_mi2hed->SetDataFolder(DataFolder);

    m_vPatientIds = m_p_mi2hed->GetPatientIDs();

    //Get reference image to set the network's dimensions
    if(!m_p_mi2hed->LoadVolumesForPatient(m_vPatientIds[0])){
        std::cerr << "Error: Could not load reference volume" << std::endl;
    }

    const ImageType3D::SizeType & clSize = m_p_mi2hed->GetReferenceVolume()->GetBufferedRegion().GetSize();

    m_channels = m_p_mi2hed->GetImageSlices(m_vPatientIds[0], 0).size();;

    int width = clSize[0];
    int height = clSize[1];
    
    Caffe::set_mode(Caffe::GPU);
    Caffe::SetDevice(iDeviceId);
            
    // Load the network
    m_p_net.reset(new Net<float>(DeployProtoFile, TEST));
    m_p_net->CopyTrainedLayersFrom(TrainedModelFile);

    //use blob to put data into the net
    caffe::Blob<float> * p_blob_input = m_p_net->input_blobs()[0];

    p_blob_input->Reshape(1, m_channels, width, height);

    m_p_net->Reshape();

    const int output_blob_index = m_p_net->output_blob_indices()[iOutputBlobIndex];
    const string & output_name = m_p_net->blob_names()[output_blob_index];

    std::cout << "Info: " << "output blob index: " << output_blob_index << "    output name: " << output_name << std::endl;

    m_strOutputName = output_name;

}

bool Detector::Detect(){

    float fAverageLoss=0;

    for (int i = 0; i < m_vPatientIds.size(); ++i){
        std::string strPatient = m_vPatientIds[i];
        std::cout << "Info: Running detection for patient " << strPatient << "." << std::endl;
        
        if (!m_p_mi2hed->LoadVolumesForPatient(strPatient)){
            std::cerr << "Error: could not load volumes for patient " << strPatient << std::endl;
            return false;
        }

        const ImageType3D::SizeType &clSize = m_p_mi2hed->GetReferenceVolume()->GetBufferedRegion().GetSize();
        const ImageType3D::SpacingType clSpacing = m_p_mi2hed->GetReferenceVolume()->GetSpacing();

        ImageType3D::Pointer p_clProbMap = ImageType3D::New();
        p_clProbMap->SetRegions(clSize);
        p_clProbMap->SetSpacing(clSpacing);
        p_clProbMap->SetOrigin(m_p_mi2hed->GetReferenceVolume()->GetOrigin());
        p_clProbMap->SetDirection(m_p_mi2hed->GetReferenceVolume()->GetDirection());

        p_clProbMap->Allocate(true);
        p_clProbMap->FillBuffer(0);
        for (int z = 0; z < clSize[2]; ++z){
            if(!DetectOne(strPatient, z)){
                std::cerr << "Error: Unable to run detection for " << strPatient << "." << std::endl;
                return false;
            }
            //Get output
            const boost::shared_ptr<caffe::Blob<float> > p_blob_output = m_p_net->blob_by_name(m_strOutputName);

            const float *  output = p_blob_output->cpu_data();

            for (int y = 0; y < clSize[1]; ++y){
                for (int x = 0; x < clSize[0]; ++x){
                    ImageType3D::IndexType clIndex = {x, y, z};
                    bool bInProstate = m_p_mi2hed->GetForground()->GetPixel(clIndex) != 0;
                    if(!bInProstate)
                        continue;
                    short value = (short)(output[y*clSize[0]+x]*1024);
                    p_clProbMap->SetPixel(clIndex, value);

                }
            }

        }

        if (m_bColor){
            std::string strFolderName = "ColorPredictions";
            nih::MakeFolder(strFolderName.c_str());
            std::string strFileName = "ColorPredictions/" + strPatient + ".mhd";
            m_p_mi2hed->SetThreshold(m_fThreshold);
            if(!m_p_mi2hed->SaveColorProbabilityMap(p_clProbMap, strFileName, strPatient))
                return false;

        }
        else{
        
            std::string strFolderName = "ProbabilityMaps";
            nih::MakeFolder(strFolderName.c_str());
            std::string strFileName = "ProbabilityMaps/" + strPatient + ".mhd";
            if(!nih::SaveVolume<short>(p_clProbMap, strFileName.c_str(), true)){
                std::cerr << "Error: Unable to save volume to file " << strFileName << "." << std::endl;
                return false;
            }

        }
    }
    return true;
}

bool Detector::DetectOne(std::string &  strPatientId, int z){
    
    std::vector<itk::Image<short, 2>::Pointer> vImages = m_p_mi2hed->GetImageSlices(strPatientId, z);

    std::vector<itk::Image<short, 2>::Pointer> vCompressed(3);

    std::string strT2WI = "T2WI";
    std::string strADC = "ADC";
    std::string strB2000 = "B2000";

    vCompressed[0] = m_p_mi2hed->CompressImage(vImages[0], strT2WI);
    vCompressed[1] = m_p_mi2hed->CompressImage(vImages[1], strADC);
    vCompressed[2] = m_p_mi2hed->CompressImage(vImages[2], strB2000);
    
    caffe::Blob<float> * p_blob_input = m_p_net->input_blobs()[0];
    
    float * input_data = p_blob_input->mutable_cpu_data();

    const ImageType::SizeType & clSize = vImages[0]->GetBufferedRegion().GetSize();
    //CHECK(reinterpret_cast<float*>(vImages[0]->GetBufferPointer()) == m_p_net->input_blobs()[0]->cpu_data())
    //    << "Input channels are not wrapping the input layer of the network.";
    int N = clSize[0]*clSize[1];

    for (int i = 0; i < vImages.size(); ++i){

        std::copy(vCompressed[2-i]->GetBufferPointer(), vCompressed[2-i]->GetBufferPointer()+N, input_data);
        //std::copy(vImages[2-i]->GetBufferPointer(), vImages[2-i]->GetBufferPointer()+N, input_data);
        input_data+=N;
    }

    //Run network
    std::vector<caffe::Blob<float>*> bottom_vec;
    float f_loss;
    m_p_net->Forward(bottom_vec, &f_loss);
    
    return true;

}

void Usage(const char * p_cArg0){
  std::cerr << "Usage: " << p_cArg0 << " -d DataFolder -g GPUDeviceId -l ListFile -p ProtoFile -m TrainedModel -b BlobIndex -t Threshold c(colorPredictions) r(run)" << std::endl;
  exit(1);
}


int main(int argc, char **argv){



	const char * const p_cArg0 = argv[0];
        std::string strProtoFile, strTrainedModel, strDataFolder, strListFile;
        bool bTest;
        int iDeviceId = 0;
        bTest=false;
        int iBlobIndex=5;
        bool bColor=false;
        float fThreshold=0.95;

	int c = 0;
	while ((c = getopt(argc, argv, "b:cd:hl:g:m:p:rt:")) != -1) {
		switch(c) {
                        case 'd':
                                strDataFolder = optarg;
                                break;
                        case 'b':
                                iBlobIndex=std::atoi(optarg);
                                break;
                        case 'c':
                                bColor=true;
                                break;
                        case 'g':
                                iDeviceId = std::atoi(optarg);
                                break;
			case 'h':
				Usage(p_cArg0); // Exits
				break;
                        case 'l':
                                strListFile = optarg;
                                break;
                        case 'p':
                                strProtoFile = optarg;
                                break;
                        case 'r':
                                bTest = true;
                                break;
                        case 't':
                                fThreshold=std::atof(optarg);
                                break;
                        case 'm':
                                strTrainedModel = optarg;
                                break;
			case '?':
			default:
				Usage(p_cArg0); // Exits
		}
	}

        if(bTest){

            Detector detector(strDataFolder, strProtoFile, strTrainedModel, strListFile, iDeviceId, iBlobIndex);
            
            detector.SetColorOutput(bColor);
            detector.SetThreshold(fThreshold);

            if(!detector.Detect()){
                std::cerr << "Error: Failed to run detection " << std::endl;
                return -1;
            }

        }

        return 0;
}
