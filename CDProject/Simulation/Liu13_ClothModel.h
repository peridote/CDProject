#include "Common/Common.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include "Simulation/TriangleModel.h"


using namespace PBD;

class Liu13_ClothModel : public TriangleModel
{
	public:
		Liu13_ClothModel();
		~Liu13_ClothModel();
		enum Springtype{st, sh, be};
		struct Edge
		{
			unsigned int m_vert[2];
			Real coeff = 1;
			Springtype type;
		};

		typedef Eigen::SparseMatrix<double> SparseMatrix;
		typedef Eigen::Triplet<Real> Triplet;
		typedef Edge Spring;

	protected:
		unsigned int ncol;
		unsigned int nrow;
		Real m_mass;
		SparseMatrix m_massMatrix;
		SparseMatrix m_L;
		SparseMatrix m_J;
		SparseMatrix m_ML;
		Eigen::VectorXd m_X;
		Eigen::VectorXd m_XOld;
		Eigen::VectorXd m_d;
		Eigen::VectorXd m_b;
		std::vector<Real> m_restLength;
		std::vector<Spring> m_springs;
		Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
		std::set<unsigned int> m_fixed_points;

		int size;

	public:
		void getPositionVector(ParticleData& pd);
		void getMassMatrix(Eigen::MatrixXd& matrix);
		void setMassMatrix(ParticleData& pd);
		void setParticleMass(ParticleData& pd);
		void setLMatrix();
		void setJMatrix();
		void setMLMatrix(Real h);
		void setRestLength(ParticleData& pd);
		void getForceVector(ParticleData& pd, Real h);
		void setVectorSize();
		void localStep(ParticleData& pd);
		void globalStep(ParticleData& pd, Real h);
		void vectorToPosition(ParticleData& pd);
		void addSprings();
		void addStructureSpring();
		void addShearSpring();
		void addBendingSpring();
		void setFixPoint(ParticleData& pd);
		void fixOverSpring(ParticleData& pd, Real h2);
		void collisionSphere(ParticleData& pd);

		void fixedPointMovement(ParticleData& pd);


		FORCE_INLINE Real getMass()
		{
			return m_mass;
		}
		std::vector<Spring>& getSprings() { return m_springs; }
};


