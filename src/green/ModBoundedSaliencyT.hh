#ifndef OPENMESH_DECIMATER_MODBOUNDEDSALIENCYT_HH
#define OPENMESH_DECIMATER_MODBOUNDEDSALIENCYT_HH


//== INCLUDES =================================================================

#include <OpenMesh/Tools/Decimater/ModBaseT.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Geometry/NormalConeT.hh>
#include "FloatHistogram.h"

using namespace OpenMesh;
using namespace OpenMesh::Decimater;

namespace lce {

	template <class MeshT>
	class ModBoundedSaliencyT: public OpenMesh::Decimater::ModBaseT<MeshT>
	{
	public:

		// Defines the types Self, Handle, Base, Mesh, and CollapseInfo
		// and the memberfunction name()
		DECIMATING_MODULE(ModBoundedSaliencyT, MeshT, BoundedSaliency)

	public:

		/** Constructor
		*  \internal
		*/
		ModBoundedSaliencyT(MeshT& _mesh)
			: Base(_mesh, true),minCollapsePercentile(0),maxCollapsePercentile(1)
		{

		}


		/// Destructor
		virtual ~ModBoundedSaliencyT()
		{
		}


	public: // inherited

		virtual void initialize()
		{
			
			green::TriMesh::VertexIter vIt, vEnd;
			float s;
			if(saliencyProperty.is_valid()){
				for (vIt = Base::mesh().vertices_begin(), vEnd = Base::mesh().vertices_end(); vIt != vEnd; ++vIt)
				{
					s = Base::mesh().property(saliencyProperty, *vIt);
					weightHistogram_.add(s, 1.0f);
				}
			}
			else
			{
				std::cerr << "No vertex weighting histogram specified, decimation module cannot operate." << std::endl;
				this->minCollapseSaliency = std::numeric_limits<double>::infinity();
			}
		}


		virtual float collapse_priority(const CollapseInfo& _ci)
		{
		
			float saliencySum = 0.0f;

		
			//saliencySum = std::floor((Base::mesh().property(saliencyProperty, _ci.v0)+ Base::mesh().property(saliencyProperty, _ci.v1))/2.0);
			saliencySum = std::min(Base::mesh().property(saliencyProperty, _ci.v0) , Base::mesh().property(saliencyProperty, _ci.v1));

			//saliencySum = Base::mesh().property(saliencyProperty, _ci.v1);
			return ((saliencySum > this->minCollapseSaliency&& saliencySum <= this->maxCollapseSaliency) ?
				float(Base::LEGAL_COLLAPSE) : float(Base::ILLEGAL_COLLAPSE));
			
		}


		virtual void postprocess_collapse(const CollapseInfo& _ci)
		{
			////TODO: optimize
			//OpenMesh::VPropHandleT<float> saliencyProperty;
			//if (Base::mesh().get_property_handle(saliencyProperty, "quality"))
			//{
			//	Base::mesh().property(saliencyProperty, _ci.v1) += Base::mesh().property(saliencyProperty, _ci.v0);
			//	//here's the magic factor! can be used to steer triangle density
			//	Base::mesh().property(saliencyProperty, _ci.v1) *= this->error_tolerance_factor_;
			//}

		}


	public: // specific methods



		void setSaliencyProperty(OpenMesh::VPropHandleT<float> prop)
		{
			this->saliencyProperty = prop;
		}

		void setMinCollapsePercentile(float min)
		{
			this->minCollapsePercentile = min;
			this->minCollapseSaliency = this->weightHistogram_.getPercentile(minCollapsePercentile);
			std::cout << "min collapse saliency is " << minCollapseSaliency << std::endl;


		}

		void setMaxCollapsePercentile(float max)
		{
			this->maxCollapsePercentile = max;
			this->maxCollapseSaliency = this->weightHistogram_.getPercentile(maxCollapsePercentile);
			std::cout << "max collapse saliency is " << maxCollapseSaliency << std::endl;

		}

	private:

		OpenMesh::VPropHandleT<float> saliencyProperty;

		lce::FloatHistogram weightHistogram_;
		float minCollapsePercentile;
		float maxCollapsePercentile;

		float minCollapseSaliency;
		float maxCollapseSaliency;


	};
}

#endif // OPENMESH_DECIMATER_MODBOUNDEDSALIENCYT_HH defined