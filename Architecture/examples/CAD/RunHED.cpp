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

class Trainer {
    public:
        Trainer(const std::string & DataFolder,
                const std::string & ProtoFile,
                const std::string & TrainedModelFile,
                const std::string & ListFile,
                int iDeviceId);

        bool Train();
        bool CreateDB();
        bool WriteDB(const int iX,
                    const int iY,
                    const int iBatchSize,
                    const int iChannels,
                    const int RANK );
        bool ReadDB(const int iX,
                    const int iY,
                    const int iBatchSize,
                    const int iChannels,
                    const int RANK );

        typedef itk::Image<short, 2> ImageType;
        typedef itk::Image<short, 3> ImageType3D;

        virtual ~Trainer() { }

    private:
        std::shared_ptr<caffe::Net<float> > m_p_net;
        nih::mi2hed * m_p_mi2hed;
        int m_channels;
        std::vector<std::string> m_vPatientIds;
        std::vector<ImageType::Pointer> m_vSlices;
};


Trainer::Trainer(const std::string & DataFolder,
                const std::string & ProtoFile,
                const std::string & TrainedModelFile,
                const std::string & ListFile,
                int iDeviceId){

    m_p_mi2hed = new nih::mi2hed;

    m_p_mi2hed->LoadPatientIds(ListFile);

    m_p_mi2hed->SetDataFolder(DataFolder);

    m_vPatientIds = m_p_mi2hed->GetPatientIDs();

    //if(!CreateDB())
      //  std::cerr << "Error: Unable to create database" << std::endl;

    //Get reference image to set the network's dimensions
    //std::cout << "Info: Loading reference volume to reshape layers" << std::endl;
    //if(!m_p_mi2hed->LoadVolumesForPatient(m_vPatientIds[0])){
    //    std::cerr << "Error: Could not load reference volume" << std::endl;
    //}

    //const ImageType3D::SizeType & clSize = m_p_mi2hed->GetReferenceVolume()->GetBufferedRegion().GetSize();

    //m_channels = m_p_mi2hed->GetImageSlices(m_vPatientIds[0], 0).size();;

    //int width = clSize[0];
    //int height = clSize[1];
    
    //Caffe::set_mode(Caffe::GPU);
    //Caffe::SetDevice(iDeviceId);
            
    // Load the network
    //m_p_net.reset(new Net<float>(ProtoFile, TRAIN));
    //m_p_net->CopyTrainedLayersFrom(TrainedModelFile);

}

bool Trainer::CreateDB(){

    std::cout << "Info: bla bla " << std::endl;

    //Load reference volume for database dimensions
    std::cout << "Info: Loading reference volume for database dimensions" << std::endl;
    if(!m_p_mi2hed->LoadVolumesForPatient(m_vPatientIds[0])){
        std::cerr << "Error: Could not load reference volume" << std::endl;
    }

    const ImageType3D::SizeType & clSize = m_p_mi2hed->GetReferenceVolume()->GetBufferedRegion().GetSize();

    const int iX = clSize[0];
    const int iY = clSize[1];
    const int iBatchSize = 1;
    const int iChannels = 3;
    const int RANK = 4;

    if(!WriteDB(iX, iY, iBatchSize, iChannels, RANK))
        return false;

    if(!ReadDB(iX, iY, iBatchSize, iChannels, RANK))
        return false;

    return true;
}

bool Trainer::WriteDB(const int iX,
            const int iY,
            const int iBatchSize,
            const int iChannels,
            const int RANK ){

    std::vector<itk::Image<short, 2>::Pointer> vImages = m_p_mi2hed->GetImageSlices(m_vPatientIds[0], 13);

    float fImage[iBatchSize][iChannels][iX][iY];
    for (int c = 0; c < iChannels; c++){
        for (int y = 0; y < iY; y++){
            for (int x = 0; x < iX; x++){
                const ImageType::IndexType & clIndex = {x, y};
                fImage[0][c][x][y] = (float)vImages[c]->GetPixel(clIndex);
            }
        }
    }
    try
    {
        Exception::dontPrint();

        const H5std_string FILE_NAME("test");
        const H5std_string DATASET_NAME("Image");
        const H5std_string DATASET_NAME2("GroundTruth");
        H5File file (FILE_NAME, H5F_ACC_TRUNC );

        hsize_t dimsf[4];
        dimsf[0] = iBatchSize;
        dimsf[1] = iChannels;
        dimsf[2] = iX;
        dimsf[3] = iY;

        DataSpace dataspace( RANK, dimsf );

        H5::FloatType datatype (PredType::NATIVE_FLOAT);
        datatype.setOrder(H5T_ORDER_LE);


        DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
        
        dataset.write( fImage, PredType::NATIVE_FLOAT);
        
    }

    catch( FileIException error )
    {
        error.printError();
        return false;
    }

    catch ( DataSetIException error)
    {
        error.printError();
        return false;
    }
    catch (DataSpaceIException error )
    {
        error.printError();
        return false;
    }
    catch (DataTypeIException error)
    {
        error.printError();
        return false;
    }

    return true;
}


bool Trainer::ReadDB(const int iX,
            const int iY,
            const int iBatchSize,
            const int iChannels,
            const int RANK ){

    const int NX_SUB = 3;
    const int NY_SUB = 4;
    const int RANK_OUT = RANK;
    //Output buffer initialization.
    float fImage[iBatchSize][iChannels][iX][iY];
    for (int c = 0; c < iChannels; c++){
        for (int y = 0; y < iY; y++){
            for (int x = 0; x < iX; x++){
                fImage[0][c][x][y] = 0;
            }
        }
    }

    try
    {
        Exception::dontPrint();

        H5File file ("test", H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet("Image");

        H5T_class_t type_class = dataset.getTypeClass();

        if( type_class == H5T_FLOAT )
            std::cout << "Data set has FLOAT type" << std::endl;

        DataSpace dataspace = dataset.getSpace();

        int rank = dataspace.getSimpleExtentNdims();

        hsize_t dims_out[2];
        int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);

        std::cout << "rank " << rank << ", dimensions " << 
            (unsigned long)(dims_out[0]) << " x " <<
            (unsigned long)(dims_out[1]) << " x " <<
            (unsigned long)(dims_out[2]) << " x" <<
            (unsigned long)(dims_out[3]) << std::endl;

        hsize_t offset[2];
        hsize_t count[2];
        offset[0] = 1;
        offset[1] = 2;
        count[0] = 1;
        count[1] = 1;

        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset );

        hsize_t dimsm[4];
        dimsm[0] = iBatchSize;
        dimsm[1] = iChannels;
        dimsm[2] = iX;
        dimsm[3] = iY;

        DataSpace memspace(RANK_OUT, dimsm );

        hsize_t offset_out[3];   // hyperslab offset in memory
        hsize_t count_out[3];    // size of the hyperslab in memory
        offset_out[0] = 3;
        offset_out[1] = 0;
        offset_out[2] = 0;
        count_out[0]  = NX_SUB;
        count_out[1]  = NY_SUB;
        count_out[2]  = 1;
        memspace.selectHyperslab( H5S_SELECT_SET, count_out, offset_out );

        dataset.read(fImage, PredType::NATIVE_FLOAT, memspace, dataspace);

        ImageType::Pointer p_clImage = ImageType::New();

        ImageType::SizeType clSize = {iX, iY};

    }

    catch ( FileIException error )
    {
        error.printError();
        return false;
    }

    catch ( DataSetIException error) 
    {
        error.printError();
        return false;
    }

    catch ( DataSpaceIException error )
    {
        error.printError();
        return false;
    }

    catch ( DataTypeIException error )
    {
        error.printError();
        return false;
    }
    return true;
}

bool Trainer::Train() {
    return CreateDB();
}



class Detector {

    public:
        Detector(const std::string & DataFolder,
                const std::string & DeployProtoFile,
                const std::string & TrainedModelFile, 
                const std::string & ListFile,
                int iDeviceId);
        
        bool Detect();
        bool DetectOne(std::string & strPatientId, int z);

        typedef itk::Image<short, 2> ImageType;
        typedef itk::Image<short, 3> ImageType3D;

        virtual ~Detector() { }

    private:
        std::shared_ptr<caffe::Net<float> > m_p_net;
        nih::mi2hed * m_p_mi2hed;
        std::vector<std::string> m_vPatientIds;
        int m_channels;
        std::string m_strOutputName;
        float m_floss;
};

Detector::Detector(const string & DataFolder,
                    const std::string & DeployProtoFile,
                    const std::string & TrainedModelFile, 
                    const std::string & ListFile,
                    int iDeviceId){
    
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

    const int output_blob_index = m_p_net->output_blob_indices()[5];
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
            fAverageLoss += m_floss;

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
        
        std::string strFolderName = "ProbabilityMaps";

        nih::MakeFolder(strFolderName.c_str());

        std::string strFileName = "ProbabilityMaps/" + strPatient + ".mhd";

        if(!nih::SaveVolume<short>(p_clProbMap, strFileName.c_str(), true)){
            std::cerr << "Error: Unable to save volume to file " << strFileName << "." << std::endl;
            return false;
        }

        std::cout << "Average Loss for Patient: " << fAverageLoss << std::endl;

        fAverageLoss = 0;
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
    
    m_floss = f_loss;

    return true;

}

void Usage(const char * p_cArg0){
  std::cerr << "Usage: " << p_cArg0 << " -d DataFolder -i inputFileName -g GPUDeviceId -l ListFile -p ProtoFile -o OutputFileName -m TrainedModel [rRS]" << std::endl;
  exit(1);
}


int main(int argc, char **argv){



	const char * const p_cArg0 = argv[0];
        std::string strProtoFile, strTrainedModel, strDataFolder, strListFile, strOutputFileName, strFileName;
        bool bTest, bTrain, bRunPNGFile;
        int iDeviceId = 0;
        bTest=bTrain=bRunPNGFile=false;

	int c = 0;
	while ((c = getopt(argc, argv, "d:hi:l:g:m:o:p:rRSt:")) != -1) {
		switch(c) {
                        case 'd':
                                strDataFolder = optarg;
                                break;
                        case 'g':
                                iDeviceId = std::atoi(optarg);
                                break;
			case 'h':
				Usage(p_cArg0); // Exits
				break;
                        case 'i':
                                strFileName = optarg;
                                break;
                        case 'l':
                                strListFile = optarg;
                                break;
                        case 'p':
                                strProtoFile = optarg;
                                break;
                        case 'o':
                                strOutputFileName = optarg;
                                break;
                        case 'r':
                                bTest = true;
                                break;
                        case 'R':
                                bTrain = true;
                                break;
                        case 'S':
                                bRunPNGFile = true;
                                break;
                        case 'm':
                                strTrainedModel = optarg;
                                break;
			case '?':
			default:
				Usage(p_cArg0); // Exits
		}
	}

        if(bTrain){

            Trainer trainer(strDataFolder, strProtoFile, strTrainedModel, strListFile, iDeviceId);

            if(!trainer.Train()){
                std::cerr << "Error: Failed to train. " << std::endl;
                return -1;
            }
        }
        if(bTest){

            Detector detector(strDataFolder, strProtoFile, strTrainedModel, strListFile, iDeviceId);

            if(!detector.Detect()){
                std::cerr << "Error: Failed to run detection " << std::endl;
                return -1;
            }

        }

        if (bRunPNGFile){
            //Set GPU
            Caffe::set_mode(Caffe::GPU);
            Caffe::SetDevice(0);
            
            typedef itk::RGBPixel<unsigned char> RGBPixelType;
            typedef itk::Image<RGBPixelType, 2> RGBImageType;
            typedef itk::Image<unsigned char, 2> ImageType;
            typedef itk::ImageFileReader<RGBImageType> ReaderType;

            ReaderType::Pointer p_clReader = ReaderType::New();
            p_clReader->SetFileName(strFileName.c_str());

            try 
            {
                p_clReader->Update();
            }
            catch (itk::ExceptionObject & e)
            {
                std::cerr << "Error: " << e << std::endl;
                return -1;
            }

            RGBImageType::Pointer p_clImage = p_clReader->GetOutput();

            ImageType::Pointer clT2WI = ImageType::New();
            ImageType::Pointer clADC = ImageType::New();
            ImageType::Pointer clB2000 = ImageType::New();

            const RGBImageType::SizeType &clRGBSize = p_clImage->GetBufferedRegion().GetSize();

            clT2WI->SetRegions(clRGBSize);
            clADC->SetRegions(clRGBSize);
            clB2000->SetRegions(clRGBSize);

            clT2WI->Allocate();
            clADC->Allocate();
            clB2000->Allocate();

            for (int y=0; y<clRGBSize[1]; y++){
                for (int x=0; x<clRGBSize[0]; x++){
                    RGBImageType::IndexType clIndex = { x, y};

                    clT2WI->SetPixel(clIndex, p_clImage->GetPixel(clIndex)[0]);
                    clADC->SetPixel(clIndex, p_clImage->GetPixel(clIndex)[1]);
                    clB2000->SetPixel(clIndex, p_clImage->GetPixel(clIndex)[2]);
                }
            }

            std::vector<ImageType::Pointer> vChannels(3);
            vChannels[2] = clT2WI;
            vChannels[1] = clADC;
            vChannels[0] = clB2000;

            // Load the network
            shared_ptr<Net<float> > p_net;
            p_net.reset(new Net<float>(strProtoFile, TEST));
            p_net->CopyTrainedLayersFrom(strTrainedModel);

            //use blob to put data into the net
            caffe::Blob<float> * p_blob_input = p_net->input_blobs()[0];
            int inum_channels = p_blob_input->channels();

            p_blob_input->Reshape(1, 3, 512,512);

            p_net->Reshape();


            //Put images into the input layer
            int width = p_blob_input->width();
            int height = p_blob_input->height();
            float * input_data = p_blob_input->mutable_cpu_data();
            for (int i=0; i < p_blob_input->channels(); ++i) {

                std::copy(vChannels[i]->GetBufferPointer(), vChannels[i]->GetBufferPointer() + width*height, input_data);
                input_data += width * height;
            }

            //Run network
            std::vector<caffe::Blob<float>*> bottom_vec;
            float f_loss;
            p_net->Forward(bottom_vec, &f_loss);


            //Get output
            const int output_blob_index = p_net->output_blob_indices()[5];
            const string & output_name = p_net->blob_names()[output_blob_index];

            std::cout << "Info: " << "output blob index: " << output_blob_index << "    output name: " << output_name << std::endl;

            const boost::shared_ptr<caffe::Blob<float> > p_blob_output = p_net->blob_by_name(output_name);

            const float *  output = p_blob_output->cpu_data();

            //Convert image output

            ImageType::Pointer p_image = ImageType::New();

            int W = p_blob_output->width();
            int H = p_blob_output->height();
            int C = p_blob_output->channels();

            std::cout << W << " " << H << " " << C << std::endl;

            ImageType::SizeType clSize;
            clSize[0] = W;
            clSize[1] = H;
            p_image->SetRegions(clSize);
            p_image->Allocate(true);

            for (int h = 0; h < H; ++h){
                for ( int w = 0; w < W; ++w){
                    ImageType::PixelType pixel;

                    ImageType::IndexType index = {w,h};

                    pixel = (unsigned char)(*(output+h*W+w)*255);
                    p_image->SetPixel(index, pixel);

                    }
            }
            
            //Save image to file
            if(!nih::SaveSlice<unsigned char>(p_image, strOutputFileName.c_str())){
                std::cerr << "Error: Failed to save file." << std::endl;
                return -1;
            }

        }

        return 0;
}

