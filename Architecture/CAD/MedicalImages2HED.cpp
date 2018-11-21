#include "MedicalImages2HED.h"
#include "Common.h"
#include "VoiFileIO.h"
#include "VoiToMask.h"
#include "itkRescaleIntensityImageFilter.h"
#include <unordered_set>
#include "itkJetColormapFunction.h"

namespace nih {

	bool mi2hed::LoadPatientIds(const std::string &strListFile){
		m_vPatientIds.clear();

		std::ifstream ListStream(strListFile.c_str());

		if (!ListStream) {
			std::cerr << "Error: Could not load '" << strListFile << "'." << std::endl;
			return false;
		}

		std::string strLine;
		while (std::getline(ListStream, strLine)) {
			size_t p = strLine.find('\r');

			if (p != std::string::npos)
				strLine.erase(p);

			if (strLine.empty())
				continue;

		m_vPatientIds.push_back(strLine);
      }

		return m_vPatientIds.size() > 0;
    }

	template<typename PixelType>
	bool mi2hed::ConvertImage(itk::Image<short, 2>::Pointer &p_clImage, typename itk::Image<PixelType, 2>::Pointer &p_clOutputImage, std::string strCodeBookFile /*T2WI | ADC | B2000*/){
	  
		/*typename*/ CodeBookTree<short,PixelType> clCodeBookTree;

		if (!/*typename*/ clCodeBookTree.LoadFromFile(strCodeBookFile.c_str())) {
			std::cerr << "Error: Could not load Codebook File. " << strCodeBookFile << std::endl;
			return false;
		}
		p_clOutputImage->SetRegions(p_clImage->GetBufferedRegion().GetSize());
		p_clOutputImage->SetSpacing(p_clImage->GetSpacing());
		p_clOutputImage->SetOrigin(p_clImage->GetOrigin());
		p_clOutputImage->SetDirection(p_clImage->GetDirection());
		p_clOutputImage->SetMetaDataDictionary(p_clImage->GetMetaDataDictionary());

		p_clOutputImage->Allocate();

		const short * const p_inBuffer = p_clImage->GetPixelContainer()->GetBufferPointer();
		const size_t length = p_clImage->GetPixelContainer()->Size();
		PixelType * const p_outBuffer = p_clOutputImage->GetPixelContainer()->GetBufferPointer();

		for (size_t i = 0; i < length; ++i) {
			if (!clCodeBookTree.Convert(p_inBuffer[i], p_outBuffer[i])) {
			  std::cerr << "Error: Conversion failed." << std::endl;
			  return false;
			}
		}

		return true;
	}

        
        itk::Image<short, 2>::Pointer mi2hed::CompressImage(itk::Image<short, 2>::Pointer &p_clImage, std::string strImageSequenceType){

		//Convert images
		std::string strCodeBook = GetDataFolder() + "/codebook" + strImageSequenceType + ".clf";

		itk::Image<short, 2>::Pointer p_clOutputImage = itk::Image<short, 2>::New();

		if (!ConvertImage<short>(p_clImage, p_clOutputImage, strCodeBook)) 
			std::cerr << "Error: Could not convert images. " << std::endl;

                return p_clOutputImage;
	
        }

	bool mi2hed::CombineImages(itk::Image<short, 2>::Pointer &p_clImage1, itk::Image<short, 2>::Pointer &p_clImage2, itk::Image<short, 2>::Pointer &p_clImage3, std::string strImageName){
		const itk::Image<short, 2>::SizeType &clSize = p_clImage1->GetBufferedRegion().GetSize();
		const itk::Image<short, 2>::SizeType &clSize2 = p_clImage2->GetBufferedRegion().GetSize();
		const itk::Image<short, 2>::SizeType &clSize3 = p_clImage3->GetBufferedRegion().GetSize();
#ifdef CODEBOOK    
		//Convert images
		std::string strCodeBook1 = GetModelFolder() + "/codebookT2WI.clf";
		std::string strCodeBook2 = GetModelFolder() + "/codebookADC.clf";
		std::string strCodeBook3 = GetModelFolder() + "/codebookB2000.clf";

		itk::Image<short, 2>::Pointer p_clOutputImage1 = itk::Image<short, 2>::New();
		itk::Image<short, 2>::Pointer p_clOutputImage2 = itk::Image<short, 2>::New();
		itk::Image<short, 2>::Pointer p_clOutputImage3 = itk::Image<short, 2>::New();

		if (!ConvertImage<short>(p_clImage1, p_clOutputImage1, strCodeBook1) || 
			!ConvertImage<short>(p_clImage2, p_clOutputImage2, strCodeBook2) || 
			!ConvertImage<short>(p_clImage3, p_clOutputImage3, strCodeBook3) )
			std::cerr << "Error: Could not convert images. " << std::endl;
#else
	
    //Scale down image to 0-255
		for (unsigned int y = 0; y < clSize[1]; ++y){
		  for (unsigned int x = 0; x < clSize[0]; ++x){
		    itk::Image<short, 2>::IndexType clIndex = {x , y};

		    float fValue = (float)p_clImage1->GetPixel(clIndex);

		        short sSetValue = (short)std::min(255.0, 255.0*fValue/1000.0);
		        p_clImage1->SetPixel(clIndex, sSetValue);

		        fValue = (float)p_clImage2->GetPixel(clIndex);
		        sSetValue = (short)std::min(255.0, 255.0*fValue/3000.0);
		        p_clImage2->SetPixel(clIndex, sSetValue);

		        fValue = (float)p_clImage3->GetPixel(clIndex);
		        sSetValue = (short)std::min(255.0, 255.0*fValue/700.0);
		        p_clImage3->SetPixel(clIndex, sSetValue);
		  }
		}
#endif	
		ComposeRGBFilterType::Pointer p_clComposeRGB = ComposeRGBFilterType::New();
		MakeFolder(GetCombinedImageFolder().c_str());
		const std::string strOutputName = GetCombinedImageFolder() + '/' + strImageName + ".png";
#ifdef CODEBOOK
		p_clComposeRGB->SetInput1(p_clOutputImage1);
		p_clComposeRGB->SetInput2(p_clOutputImage2);
		p_clComposeRGB->SetInput3(p_clOutputImage3);
#else
		p_clComposeRGB->SetInput1(p_clImage1);
		p_clComposeRGB->SetInput2(p_clImage2);
		p_clComposeRGB->SetInput3(p_clImage3);
#endif
		try
		{
		  p_clComposeRGB->Update();
		}
		catch(itk::ExceptionObject &e)
		{
		  std::cerr << e.GetDescription() << std::endl;
		  return NULL;
		}

		RGBType::Pointer p_clRGBImage = p_clComposeRGB->GetOutput();

		itk::ImageFileWriter<RGBType>::Pointer writer = itk::ImageFileWriter<RGBType>::New();
		writer->SetFileName(strOutputName);
		writer->SetInput(p_clRGBImage);
		try
		{
		  writer->Update();
		}
		catch(itk::ExceptionObject &e)
		{
		  std::cerr << e.GetDescription() << std::endl;
		}

		return true;
	}

        std::vector<itk::Image<short, 2>::Pointer> mi2hed::GetImageSlices(std::string &strPatientId, int z){

		//Check if images are loaded
		if(!m_p_clT2WIVolume){
		  std::cerr << "Error: Volumes have not been loaded for patient " << strPatientId << "." << std::endl;
		}

		SliceContext clSlice;
		GetSliceContext(clSlice, z);

                std::vector<itk::Image<short, 2>::Pointer> vImages;

                vImages.push_back(clSlice.p_clT2WISlice);
                vImages.push_back(clSlice.p_clADCSlice);
                vImages.push_back(clSlice.p_clB2000Slice);

                return vImages;

        }

	bool mi2hed::Run() {

		MakeFolder(GetOutputFolder().c_str());

		for (int i = 0; i<m_vPatientIds.size(); ++i){
			std::string strPatient = m_vPatientIds[i];
			std::cout << "Inferring for Patient: " << strPatient << std::endl;

			//Load Images
			if(!LoadVolumesForPatient(strPatient)){
			  std::cerr << "Error: Could not load volumes for patient " << strPatient << "." << std::endl;
			  return false;
			}

			const ImageType::SizeType &clSize = m_p_clT2WIVolume->GetBufferedRegion().GetSize();
			const ImageType::SpacingType &clSpacing = m_p_clT2WIVolume->GetSpacing();

			int W, H;

			W = clSize[0];
			H = clSize[1];

			for ( int z = 0; z<clSize[2]; ++z){
				//Get slices
				SliceContext clSlice;
				GetSliceContext(clSlice, z);
				std::string strImageName=strPatient + "." + std::to_string((int)z);
				std::string GTFolder = GetGTImageFolder() + '/';
				MakeFolder(GTFolder.c_str());
				std::string strOutputName = GTFolder + strImageName + ".png";

				RGBType::Pointer p_clRGBGroundT = RGBType::New();
				p_clRGBGroundT->SetRegions(clSlice.p_clTumorSegSlice->GetBufferedRegion().GetSize());
				p_clRGBGroundT->SetSpacing(clSlice.p_clTumorSegSlice->GetSpacing());
				p_clRGBGroundT->SetOrigin(clSlice.p_clTumorSegSlice->GetOrigin());
				p_clRGBGroundT->SetDirection(clSlice.p_clTumorSegSlice->GetDirection());

				p_clRGBGroundT->Allocate(true);
				
				for (int y = 0; y < clSize[1]; ++y){
					for (int x = 0; x < clSize[0]; ++x){
						itk::RGBPixel<unsigned char> clPixel;
						RGBType::IndexType clIndex = { x, y};
						bool bTumor = clSlice.p_clTumorSegSlice->GetPixel(clIndex) != 0;
                                                bool bProstate = clSlice.p_clProstateSlice->GetPixel(clIndex) != 0;

						if(bTumor){
							clPixel[0]=255;
							clPixel[1]=255;
							clPixel[2]=255;
							p_clRGBGroundT->SetPixel(clIndex, clPixel);
						}
						else {
							clPixel[0]=0;
							clPixel[1]=0;
							clPixel[2]=0;
							p_clRGBGroundT->SetPixel(clIndex, clPixel);
						}
#if 1
                                                if (!bProstate){
                                                    clSlice.p_clT2WISlice->SetPixel(clIndex, 0);
                                                    clSlice.p_clADCSlice->SetPixel(clIndex, 0);
                                                    clSlice.p_clB2000Slice->SetPixel(clIndex, 0);
#endif //Crop end
                                                }
					}
				}
                                // Generate Combined Images
				CombineImages(clSlice.p_clT2WISlice, clSlice.p_clADCSlice, clSlice.p_clB2000Slice, strImageName);

                                // Generate GroundTruth Images
				itk::ImageFileWriter<RGBType>::Pointer writer = itk::ImageFileWriter<RGBType>::New();
				writer->SetFileName(strOutputName);
				writer->SetInput(p_clRGBGroundT);
				try
				{
				  writer->Update();
				}
				catch(itk::ExceptionObject &e)
				{
				  std::cerr << e.GetDescription() << std::endl;
				}

			}

		}
	}

	bool mi2hed::LoadVolumesForPatient(std::string &strPatientId) {
		std::string strT2WIFilePath = GetT2WIDataFolder() + '/' + strPatientId + ".%d.dcm";
		std::string strADCFilePath = GetADCDataFolder() + '/' + strPatientId + ".%d.dcm";
		std::string strB2000FilePath = GetB2000DataFolder() + '/' + strPatientId + ".%d.dcm";
		std::string strTumorMaskFolderPath = GetTumorFolder() + '/' + strPatientId;

		m_p_clT2WIVolume = LoadDicomVolumeByIndex(strT2WIFilePath.c_str());
		m_p_clADCVolume = LoadDicomVolumeByIndex(strADCFilePath.c_str());
		m_p_clB2000Volume = LoadDicomVolumeByIndex(strB2000FilePath.c_str());
		m_p_clTumorSegVolume = LoadVOITumors<unsigned char>(strTumorMaskFolderPath.c_str());

                //Load Prostate
                std::vector<unsigned char> vProstateSlices;
                m_p_clProstateVolume = itk::Image<unsigned char, 3>::New();
                m_p_clProstateVolume->SetRegions(m_p_clT2WIVolume->GetBufferedRegion().GetSize());
                m_p_clProstateVolume->SetOrigin(m_p_clT2WIVolume->GetOrigin());
                m_p_clProstateVolume->SetSpacing(m_p_clT2WIVolume->GetSpacing());
                m_p_clProstateVolume->Allocate(true);

                if(!LoadProstate(strPatientId, m_p_clProstateVolume, vProstateSlices)) {
                    std::cerr << "Error: Couldn't load prostate for patient " << strPatientId << std::endl;
                    return false;
                }
                //End of Load Prostate
                
		if(!m_p_clT2WIVolume || !m_p_clADCVolume || !m_p_clB2000Volume)
			return false;

		return true;
	}

	template<typename Type>
	typename itk::Image<Type, 3>::Pointer mi2hed::LoadVOITumors(const std::string &strFolderPath){
	  typedef itk::Image<Type, 3> TemplatedMaskType;

	  std::unordered_set<unsigned int> sSlices;

	  if (!m_p_clT2WIVolume)
	    return typename TemplatedMaskType::Pointer();

	  const ImageType::SizeType &clSize = m_p_clT2WIVolume->GetBufferedRegion().GetSize();

	  if (clSize[0] == 0 || clSize[1] == 0 || clSize[2] == 0)
	    return typename TemplatedMaskType::Pointer();

	  typename TemplatedMaskType::Pointer p_clMask = TemplatedMaskType::New();

	  std::vector<std::string> vVOIFiles;
	  FindFiles(vVOIFiles, strFolderPath.c_str(), "*.voi", false);

	  if (vVOIFiles.empty())
	    return NULL;
	  for (size_t i = 0; i < vVOIFiles.size(); ++i) {
	    const std::string &strFileName = vVOIFiles[i];

	    VoiFile clVOI;

	    if (!clVOI.LoadFromFile(strFileName.c_str()))
	      return typename TemplatedMaskType::Pointer();

	    if (!p_clMask)
	      return typename TemplatedMaskType::Pointer();

	    p_clMask->SetRegions(clSize);
	    p_clMask->SetSpacing(m_p_clT2WIVolume->GetSpacing());

	    p_clMask->Allocate(true);

	    const std::vector<VoiFile::Slice> &vVOISlices = clVOI.GetSlices();

	    for (size_t j = 0; j < vVOISlices.size(); ++j) {
	      const VoiFile::Slice &stSlice = vVOISlices[j];
	      sSlices.insert(stSlice.uiSliceNumber);
	    }

	    if (!ConvertVoiToMask3D<Type>(p_clMask, clVOI, false))
	      return typename TemplatedMaskType::Pointer();
	  }
	  return p_clMask;
	}

	void mi2hed::GetSliceContext(SliceContext &clSliceContext, unsigned int z){
	  clSliceContext.bHasT2WISlice = GetSlice<short>(clSliceContext.p_clT2WISlice, m_p_clT2WIVolume, z);
	  clSliceContext.bHasADCSlice = GetSlice<short>(clSliceContext.p_clADCSlice, m_p_clADCVolume, z);
	  clSliceContext.bHasB2000Slice = GetSlice<short>(clSliceContext.p_clB2000Slice, m_p_clB2000Volume, z);
	  clSliceContext.bHasTumorSegSlice = GetSlice<unsigned char>(clSliceContext.p_clTumorSegSlice, m_p_clTumorSegVolume, z);
          clSliceContext.bHasProstateSlice = GetSlice<unsigned char>(clSliceContext.p_clProstateSlice, m_p_clProstateVolume, z);
	}


        bool mi2hed::GenerateMHDs(){

            std::cout << "Info: Generating volumes" << std::endl;     
            
            for (int i = 0; i < GetPatientIDs().size(); ++i){

                std::string strPatientId = GetOutputFolder() + "/" + GetPatientIDs()[i] + ".%d.png";

                //Load reference volumes
                LoadVolumesForPatient(GetPatientIDs()[i]);
                const ImageType::SizeType &clSize = m_p_clT2WIVolume->GetBufferedRegion().GetSize();
                const ImageType::SpacingType &clSpacing = m_p_clT2WIVolume->GetSpacing();

                std::string strHEDOutputFolder = GetOutputFolder();

                typedef itk::NumericSeriesFileNames NamesGeneratorType;
                typedef itk::Image<short, 3> ImageType3D;
                typedef itk::ImageSeriesReader<ImageType3D> SeriesReaderType;
                typedef itk::ImageFileWriter<ImageType3D> WriterType;

                NamesGeneratorType::Pointer p_clNamesGenerator = NamesGeneratorType::New();
                p_clNamesGenerator->SetSeriesFormat( strPatientId.c_str() );
                p_clNamesGenerator->SetStartIndex( 0 );
                p_clNamesGenerator->SetEndIndex( clSize[2]-1 );
                p_clNamesGenerator->SetIncrementIndex( 1 );

                SeriesReaderType::Pointer p_clImageSeriersReader = SeriesReaderType::New();

                p_clImageSeriersReader->SetFileNames(p_clNamesGenerator->GetFileNames() );

                try
                {
                    p_clImageSeriersReader->Update();
                }
                catch (itk::ExceptionObject & e)
                {
                    std::cout << e << std::endl;
                    continue;
                }
                
                ImageType3D::Pointer p_clVolume = p_clImageSeriersReader->GetOutput();
                p_clVolume->SetOrigin(m_p_clT2WIVolume->GetOrigin());
                p_clVolume->SetSpacing(clSpacing);

#if 1 
#pragma omp parallel 
                {
#pragma omp for schedule(dynamic)
                //Scale probability map up from 256 - 1024
                    for (int z = 0; z < clSize[2]; z++){
                        for ( int y = 0; y < clSize[1]; y++){
                            for (int x = 0; x < clSize[0]; x++){
                                ImageType3D::IndexType clIndex = {x, y, z};
                                bool bProstate = m_p_clProstateVolume->GetPixel(clIndex) != 0;

                                if (!bProstate)
                                    p_clVolume->SetPixel(clIndex, 0);

                                else{

                                    short sValue = (p_clVolume->GetPixel(clIndex)*4);
                                    p_clVolume->SetPixel(clIndex, sValue);
                                }

                            }
                        }
                    }
                } //end parallel

#else
                typedef itk::RescaleIntensityImageFilter<ImageType3D, ImageType3D> RescaleType;
                RescaleType::Pointer p_clRescale = RescaleType::New();
                p_clRescale->SetInput(p_clVolume);
                p_clRescale->SetOutputMinimum(0);
                p_clRescale->SetOutputMaximum(1023);
                p_clVolume = p_clRescale->GetOutput();
#endif
                WriterType::Pointer p_clWriter = WriterType::New();
                itksys::SystemTools::MakeDirectory(m_strProbMapOutputFolder.c_str());
                std::string strOutputFile = m_strProbMapOutputFolder + "/" + GetPatientIDs()[i] + ".mhd";

                p_clWriter->SetFileName(strOutputFile.c_str() );
                                                                                                     
                p_clWriter->SetInput(p_clVolume);

                try
                {
                    p_clWriter->Update();
                }
                catch (itk::ExceptionObject & e)
                {
                    std::cout << e << std::endl;
                    continue;
                }
            }
            return true;

        }

        bool mi2hed::LoadProstate(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned char> &vSlices){
            vSlices.clear();

            std::unordered_set<unsigned char> sSlices;

            const std::string strPatientProstateFile = GetProstateFolder() + '/' + strPatient + ".voi";

            VoiFile clVOIFile;

            if (!clVOIFile.LoadFromFile(strPatientProstateFile.c_str())) {

                std::cerr << "Error: Could not load patient prostate file '" << strPatientProstateFile << "'." << std::endl;
                return false;
            }

            const std::vector<VoiFile::Slice> &vVOISlices = clVOIFile.GetSlices();

            for (size_t i = 0; i < vSlices.size(); ++i){
                const VoiFile::Slice &stSlice = vVOISlices[i];
                sSlices.insert(stSlice.uiSliceNumber);
            }

            if (!ConvertVoiToMask3D<unsigned char>(p_clMask, clVOIFile, false)) {
                std::cerr << "Error: Could not convert VOI to mask for the file '" << strPatientProstateFile << "'." << std::endl;
                return false;
            }

            vSlices.insert(vSlices.end(), sSlices.begin(), sSlices.end());
            std::sort(vSlices.begin(), vSlices.end());

            return true;

          }
    	void mi2hed::OverlayStructures(itk::Image<itk::RGBPixel<unsigned char>, 3>::Pointer p_clProbMap) {
	    if (!m_p_clProstateVolume || !m_p_clCentralGlandVolume)
		return;

	  itk::RGBPixel<unsigned char> clProstateColor, clCentralGlandColor;

	  clProstateColor.Set(126, 0, 0);
          clCentralGlandColor.Set(126, 0, 0);

	  const ImageType::SizeType &clSize = m_p_clProstateVolume->GetBufferedRegion().GetSize();
	  for (int z = 0; z < clSize[2]; ++z) {
	      for (int y = 0; y < clSize[1]; ++y) {
		  for (int x = 0; x < clSize[0]; ++x) {
			const ImageType::IndexType clIndex = { x, y, z };

			const int xBegin = std::max(0, x - 1);
			const int xEnd = std::min((int)clSize[0]-1, x+1);
			const int yBegin = std::max(0, y - 1);
			const int yEnd = std::min((int)clSize[1]-1, y+1);

			bool bProstateBoundary = false;
                        bool bCentralGlandBoundary = false;

			if (m_p_clCentralGlandVolume->GetPixel(clIndex) != 0) {
			   for (int j = yBegin; j <= yEnd && !bCentralGlandBoundary; ++j) {
				for (int i = xBegin; i <= xEnd && !bCentralGlandBoundary; ++i) {
				  const ImageType::IndexType clIndex2 = { i, j, z };

				  if (m_p_clCentralGlandVolume->GetPixel(clIndex2) <= 0)
					bCentralGlandBoundary = true;
				}
			   }
			}

			if (m_p_clProstateVolume->GetPixel(clIndex) != 0) {
			   for (int j = yBegin; j <= yEnd && !bProstateBoundary; ++j) {
				for (int i = xBegin; i <= xEnd && !bProstateBoundary; ++i) {
				  const ImageType::IndexType clIndex2 = { i, j, z };

				  if (m_p_clProstateVolume->GetPixel(clIndex2) <= 0)
					bProstateBoundary = true;
				}
			   }
			}

			if (bProstateBoundary)
			  p_clProbMap->SetPixel(clIndex, clProstateColor);
                        else if (bCentralGlandBoundary)
                            p_clProbMap->SetPixel(clIndex, clCentralGlandColor);
		  
		  }
	      }
	   }
	}

	bool mi2hed::SaveColorProbabilityMap(ImageType::Pointer &p_clVolume, std::string &strFileName, std::string &strPatientId) {
	    typedef itk::RGBPixel<unsigned char> RGBPixelType;
	    typedef itk::Image<RGBPixelType, 3> RGBImageType;
	    typedef itk::Function::JetColormapFunction<short, RGBPixelType> ColormapType;

	       if (!p_clVolume)
		   return false;
               
               const ImageType::SizeType &clSize = p_clVolume->GetBufferedRegion().GetSize();
               const ImageType::SpacingType &clSpacing = p_clVolume->GetSpacing();

               //Load Forground Prostate and Central Gland Volume
                std::vector<unsigned char> vProstateSlices;
                m_p_clProstateVolume = itk::Image<unsigned char, 3>::New();
                m_p_clProstateVolume->SetRegions(clSize);
                m_p_clProstateVolume->SetOrigin(p_clVolume->GetOrigin());
                m_p_clProstateVolume->SetSpacing(clSpacing);
                m_p_clProstateVolume->Allocate(true);

                if(!LoadProstate(strPatientId, m_p_clProstateVolume, vProstateSlices)) {
                    std::cerr << "Error: Couldn't load prostate for patient " << strPatientId << std::endl;
                    return false;
                }
                const std::string strCentralGlandPath = GetCentralGlandFolder() + "/voi/" + strPatientId + ".voi";
                m_p_clCentralGlandVolume = LoadVOIFile<unsigned char>(strCentralGlandPath.c_str());

                if(!m_p_clCentralGlandVolume)
                    std::cout << "Unable to load Central Gland Volume for patient " << strPatientId << std::endl;


	       RGBImageType::Pointer p_clRGBProbabilityMap = RGBImageType::New();

	       p_clRGBProbabilityMap->SetRegions(clSize);
	       p_clRGBProbabilityMap->SetSpacing(clSpacing);
	       p_clRGBProbabilityMap->SetOrigin(p_clVolume->GetOrigin());
	       p_clRGBProbabilityMap->SetDirection(p_clVolume->GetDirection());

	       p_clRGBProbabilityMap->Allocate(true);
               p_clRGBProbabilityMap->FillBuffer(0.0);

	       ColormapType::Pointer p_clColormap = ColormapType::New();

	       p_clColormap->SetMaximumInputValue(1024);
	       p_clColormap->SetMinimumInputValue(0);
#pragma omp parallel 
               {
#pragma omp for schedule(dynamic)
                   //Scale probability map up from 256 - 1024
                   for (int z = 0; z < clSize[2]; z++){
                       for ( int y = 0; y < clSize[1]; y++){
                           for (int x = 0; x < clSize[0]; x++){
                               ImageType::IndexType clIndex = {x, y, z};
		               if (m_p_clProstateVolume.IsNotNull() && m_p_clProstateVolume->GetPixel(clIndex) == 0)
		    	           continue;

		               short sProb = p_clVolume->GetPixel(clIndex);
                               if (sProb < m_fThreshold*1024)
                                   sProb = 0;

		               p_clRGBProbabilityMap->SetPixel(clIndex, p_clColormap->operator()(sProb));

                          }
                     }
                  }

               } //end parallel


	   OverlayStructures(p_clRGBProbabilityMap);
	   //DrawNMSPoints(p_clRGBProbabilityMap);
           //DrawGrid(p_clRGBProbabilityMap);
	   itksys::SystemTools::MakeDirectory(m_strProbMapOutputFolder.c_str());

	   SaveVolume<RGBPixelType>(p_clRGBProbabilityMap, strFileName.c_str());

        return true;
    }
template<typename DataType>
    typename itk::Image<DataType, 3>::Pointer mi2hed::LoadVOIFile(const std::string &strFilePath){
        typedef itk::Image<DataType, 3> TemplatedMaskType;

        if (!m_p_clProstateVolume)
            return typename TemplatedMaskType::Pointer();

        const ImageType::SizeType &clSize = m_p_clProstateVolume->GetBufferedRegion().GetSize();

        VoiFile clVOI;

        if(!clVOI.LoadFromFile(strFilePath.c_str()))
                return typename TemplatedMaskType::Pointer();

        typename TemplatedMaskType::Pointer p_clMask = TemplatedMaskType::New();

        if(!p_clMask)
            return typename TemplatedMaskType::Pointer();

        p_clMask->SetRegions(clSize);
        p_clMask->SetSpacing(m_p_clProstateVolume->GetSpacing());
        p_clMask->Allocate(true);

        if(!ConvertVoiToMask3D<DataType>(p_clMask, clVOI, false))
            return typename TemplatedMaskType::Pointer();

        return p_clMask;
        
    }

}
