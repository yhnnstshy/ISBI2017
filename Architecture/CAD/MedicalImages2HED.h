#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include "itkComposeImageFilter.h"
#include "CodeBookTree.h"
#include "itkRGBPixel.h"


namespace nih {

	class mi2hed {
	public:
		typedef itk::Image<short, 3> ImageType;
		typedef itk::Image<unsigned char, 3> ImageTypeUC;
		typedef itk::Image<short, 3> ProbabilityMapType;
		typedef itk::Image<unsigned char, 3> MaskType;
		typedef itk::Image<unsigned char, 2> ImageType2D;
		typedef itk::Image<itk::RGBPixel<unsigned char>, 2> RGBType;
		typedef itk::ComposeImageFilter< itk::Image<short, 2>, RGBType > ComposeRGBFilterType;
		virtual ~mi2hed() { }

		mi2hed () {
			SetDataFolder(".");
			SetOutputFolder(".");
		}

		bool LoadPatientIds(const std::string &strListFile);

		std::vector<std::string> GetPatientIDs() const {
			return m_vPatientIds;
		}

		void SetDataFolder(const std::string &strDataFolder) {
			m_strDataFolder = strDataFolder;
		}

		void SetOutputFolder(const std::string &strOutputFolder) {
			m_strOutputFolder = strOutputFolder;
		}

                void SetProbMapOutputFolder(const std::string &strProbMapOutputFolder){
                        m_strProbMapOutputFolder = strProbMapOutputFolder;
                }

                void SetThreshold(const float fThreshold){
                    m_fThreshold = fThreshold;
                }

		std::string GetDataFolder() const {
			return m_strDataFolder;
		}

		std::string GetOutputFolder() const {
			return m_strOutputFolder;
		}

		std::string GetCombinedImageFolder() const {
			return GetOutputFolder() + "/CombinedImages";
		}

		std::string GetGTImageFolder() const {
			return GetOutputFolder() + "/GroundTruth";
		}

		std::string GetT2WIDataFolder() const {
			return GetDataFolder() + "/T2WI";
		}

		// sT2Map volume data folder
		std::string GetsT2MapDataFolder() const {
			return GetDataFolder() + "/sT2Map";
		}


		// ADC volume data folder
		std::string GetADCDataFolder() const {
			return GetDataFolder() + "/ADC";
		}

		// B2000 volume data folder
		std::string GetB2000DataFolder() const {
			return GetDataFolder() + "/B2000";
		}

		std::string GetTumorFolder() const {
			return GetDataFolder() + "/Tumors";
		}


		std::string GetProstateFolder() const {
			return GetDataFolder() + "/Prostates";
		}

                ImageType::Pointer GetReferenceVolume() const {
                    return m_p_clT2WIVolume;
                }

                MaskType::Pointer GetForground() const {
                    return m_p_clProstateVolume;
                }

                std::string GetCentralGlandFolder() const {
                    return GetDataFolder() + "/CentralGlands";
                }

		bool LoadVOIs(const std::string &strFolder, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned int> &vSlices);

		template<typename Type>
		typename itk::Image<Type, 3>::Pointer LoadVOITumors(const std::string &strFilePath);

                template<typename DataType>
                typename itk::Image<DataType, 3>::Pointer LoadVOIFile(const std::string &strFilePath);

                ImageTypeUC::Pointer CreateGroundTruthFromBiopsy(std::string &strPatient);

		// This function loads all tumors for a given patient (MRN.AccessionNumber) and returns a 3D binary mask for the tumors as well as the slices where tumors exist
		bool LoadTumors(const std::string &strPatient, MaskType::Pointer p_clMask, std::vector<unsigned int> &vSlices) {
			return LoadVOIs(GetTumorFolder() + '/' + strPatient, p_clMask, vSlices);
		}

                // This function load porstate masks
                bool LoadProstate(const std::string &strPatient, itk::Image<unsigned char, 3>::Pointer p_clMask, std::vector<unsigned char> &vSlices);
		//run the whole thing
		bool Run();
		//combines three images into one RGB image
		bool CombineImages(itk::Image<short, 2>::Pointer &p_clImage1, itk::Image<short, 2>::Pointer &p_clImage2, itk::Image<short, 2>::Pointer &p_clImage3, std::string strImageName);
		//Convert images using codebook
		template <typename PixelType>
		bool ConvertImage(itk::Image<short, 2>::Pointer &p_clImage, typename itk::Image<PixelType, 2>::Pointer &p_clOutputImage, std::string strCodeBookFile);

		bool LoadVolumesForPatient(std::string &strPatientId);

                bool GenerateMHDs();

		void OverlayStructures(itk::Image<itk::RGBPixel<unsigned char>, 3>::Pointer p_clProbMap);

		bool SaveColorProbabilityMap(ImageType::Pointer &p_clVolume, std::string &strFileName, std::string &strPatientId);

                std::vector<itk::Image<short, 2>::Pointer> GetImageSlices(std::string &strPatinetId, int z);

                itk::Image<short, 2>::Pointer CompressImage(itk::Image<short, 2>::Pointer &p_clImage, std::string strImageSequenceType);

	protected:
		std::vector<std::string> m_vPatientIds;
		std::string m_strDataFolder, m_strOutputFolder, m_strProbMapOutputFolder;
                float m_fThreshold;

	private:
		struct SliceContext {
			itk::Image<short, 2>::Pointer p_clT2WISlice, p_clADCSlice, p_clB2000Slice;
			itk::Image<unsigned char, 2>::Pointer p_clTumorSegSlice, p_clProstateSlice;
			bool bHasT2WISlice, bHasADCSlice, bHasB2000Slice, bHasTumorSegSlice, bHasProstateSlice;
			SliceContext() {
				p_clT2WISlice = itk::Image<short, 2>::New();
				p_clADCSlice = itk::Image<short, 2>::New();
				p_clB2000Slice = itk::Image<short, 2>::New();
				p_clTumorSegSlice = itk::Image<unsigned char, 2>::New();
                                p_clProstateSlice = itk::Image<unsigned char, 2>::New();

				bHasT2WISlice = bHasADCSlice = bHasB2000Slice = bHasTumorSegSlice = bHasProstateSlice = false;
			}
		};
		ImageType::Pointer m_p_clT2WIVolume, m_p_clADCVolume, m_p_clB2000Volume;
		ImageTypeUC::Pointer m_p_clTumorSegVolume, m_p_clProstateVolume, m_p_clCentralGlandVolume;

		void GetSliceContext(SliceContext &clSliceContext, unsigned int z);

	};


}
