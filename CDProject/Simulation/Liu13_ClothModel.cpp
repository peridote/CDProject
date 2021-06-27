#include "Liu13_ClothModel.h"

using namespace PBD;

Liu13_ClothModel::Liu13_ClothModel()
{
	m_mass = 1.0f;
	nrow = 50;
	ncol = 50;
}

Liu13_ClothModel::~Liu13_ClothModel()
{

}


void Liu13_ClothModel::setParticleMass(ParticleData& pd)
{
	size = pd.getNumberOfParticles();
	std::vector<ParticleMesh::Face> face = m_particleMesh.getFaceData();
	std::vector<ParticleMesh::Edge> edge = m_particleMesh.getEdges();
	unsigned int vn = m_particleMesh.getVerticesPerFace();
	unsigned int numface = m_particleMesh.numFaces();
	Real pointmass = m_mass / numface;
	pointmass = pointmass / vn;
	
	for (unsigned int i = 0; i < size; i++)
	{
		pd.setMass(i, 1.0);
	}
	
	/*
	for (unsigned int i = 0; i < face.size(); i++)
	{
		for (unsigned int j = 0; j < vn; j++)
		{
			pd.setMass(edge[face[i].m_edges[j]].m_vert[0], pointmass / 2);
			pd.setMass(edge[face[i].m_edges[j]].m_vert[1], pointmass / 2);
		}
	}
	*/
}

void Liu13_ClothModel::addSprings()
{
	addStructureSpring();
	addShearSpring();
	addBendingSpring();
}

void Liu13_ClothModel::addStructureSpring()
{
	Real k = 1500;
	m_springs.reserve(m_springs.size() + 2 * (nrow-1) * (ncol-1) + nrow - 1 + ncol - 1);
	for (unsigned int i = 0; i < nrow-1; i++)
	{
		for (unsigned int j = 0; j < ncol-1; j++)
		{
			Spring spring1, spring2;
			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = i * ncol + j + 1;
			spring1.coeff = k;
			spring1.type = Springtype::st;
			m_springs.push_back(spring1);

			spring2.m_vert[0] = i * ncol + j;
			spring2.m_vert[1] = (i + 1) * ncol + j;
			spring2.coeff = k;
			spring2.type = Springtype::st;
			m_springs.push_back(spring2);
		}
	}
	for (unsigned int i = 0; i < nrow-1; i++)
	{
		Spring spring1;
		spring1.m_vert[0] = i * ncol + ncol - 1;
		spring1.m_vert[1] = (i + 1) * ncol + ncol - 1;
		spring1.coeff = k;
		spring1.type = Springtype::st;
		m_springs.push_back(spring1);
	}

	for (unsigned int i = 0; i < ncol - 1; i++)
	{
		Spring spring1;
		spring1.m_vert[0] = (nrow - 1) * ncol + i;
		spring1.m_vert[1] = (nrow - 1) * ncol + i + 1;
		spring1.coeff = k;
		spring1.type = Springtype::st;
		m_springs.push_back(spring1);
	}

}

void Liu13_ClothModel::addShearSpring()
{
	Real k = 500;
	m_springs.reserve(m_springs.size() + 2 * (nrow-1) * (ncol-1));
	for (unsigned int i = 0; i < nrow - 1; i++)
	{
		for (unsigned int j = 0; j < ncol - 1; j++)
		{
			Spring spring1, spring2;
			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = (i + 1) * ncol + j + 1;
			spring1.coeff = k;
			spring1.type = Springtype::sh;
			m_springs.push_back(spring1);

			spring2.m_vert[0] = (i + 1) * ncol + j;
			spring2.m_vert[1] = i * ncol + j + 1;
			spring2.coeff = k;
			spring2.type = Springtype::sh;
			m_springs.push_back(spring2);
		}
	}

}

void Liu13_ClothModel::addBendingSpring()
{
	Real k = 1000;
	m_springs.reserve(m_springs.size() + nrow * (ncol - 2) + (nrow - 2) * ncol);
	for (unsigned int i = 0; i < nrow; i++)
	{
		for (unsigned int j = 0; j < ncol - 2; j++)
		{
			Spring spring1;
			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = i * ncol + j + 2;
			spring1.coeff = k;
			spring1.type = Springtype::be;
			m_springs.push_back(spring1);
		}
	}
	for (unsigned int i = 0; i < nrow - 2; i++)
	{
		for (unsigned int j = 0; j < ncol; j++)
		{
			Spring spring1;

			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = (i + 2) * ncol + j;
			spring1.coeff = k;
			spring1.type = Springtype::be;
			m_springs.push_back(spring1);
		}

	}

}

void Liu13_ClothModel::setMassMatrix(ParticleData& pd)
{
	long long m;
	Real pmass;
	m = size;
	std::vector<Triplet> tripletList;
	m_massMatrix.resize(3 * static_cast<long long>(m), 3 * static_cast<long long>(m));
	for (long long i = 0; i < m; i++)
	{
		pmass = pd.getMass(i);
		tripletList.push_back(Triplet(3 * i, 3 * i, pmass));
		tripletList.push_back(Triplet(3 * i + 1, 3 * i + 1, pmass));
		tripletList.push_back(Triplet(3 * i + 2, 3 * i + 2, pmass));
	}
	m_massMatrix.setFromTriplets(tripletList.begin(), tripletList.end());
	
}

void Liu13_ClothModel::getMassMatrix(Eigen::MatrixXd& matrix)
{
	matrix = m_massMatrix;
}

void Liu13_ClothModel::init(ParticleData& pd)
{
	m_pd = &pd;
	isLiu = true;
	addSprings();
	setParticleMass(pd);

	setMassMatrix(pd);
	setLMatrix();
	setJMatrix();
	setRestLength(pd);
	setVectorSize();
	getPositionVector(pd);
	m_XOld = m_X;
	setMLMatrix(0.005);
}

void Liu13_ClothModel::getPositionVector(ParticleData& pd)
{
	long long m;
	m = size;
	m_X.setZero();

	for (long long i = 0; i < m; i++)
	{
		Vector3r pos = pd.getPosition(i + m_indexOffset);

		m_X(3 * i) = pos.x();
		m_X(3 * i + 1) = pos.y();
		m_X(3 * i + 2) = pos.z();
	}

}

void Liu13_ClothModel::setLMatrix()
{
	Real k = 1;
	std::vector<Spring> edge = m_springs;
	std::vector<Triplet> tripletList;
	m_L.resize(size*3, size*3);

	for (unsigned int i = 0; i < edge.size(); i++)
	{
		unsigned int i1 = edge[i].m_vert[0];
		unsigned int i2 = edge[i].m_vert[1];
		k = edge[i].coeff;
		//A.setZero();
		//A(i1) = 1;
		//A(i2) = -1;
		//A *= k;
		tripletList.push_back(Triplet(3 * i1, 3 * i1, k));
		tripletList.push_back(Triplet(3 * i1 + 1, 3 * i1 + 1, k));
		tripletList.push_back(Triplet(3 * i1 + 2, 3 * i1 + 2, k));

		tripletList.push_back(Triplet(3 * i1, 3 * i2, -k));
		tripletList.push_back(Triplet(3 * i1 + 1, 3 * i2 + 1, -k));
		tripletList.push_back(Triplet(3 * i1 + 2, 3 * i2 + 2, -k));

		tripletList.push_back(Triplet(3 * i2, 3 * i1, -k));
		tripletList.push_back(Triplet(3 * i2 + 1, 3 * i1 + 1, -k));
		tripletList.push_back(Triplet(3 * i2 + 2, 3 * i1 + 2, -k));

		tripletList.push_back(Triplet(3 * i2, 3 * i2, k));
		tripletList.push_back(Triplet(3 * i2 + 1, 3 * i2 + 1, k));
		tripletList.push_back(Triplet(3 * i2 + 2, 3 * i2 + 2, k));
		//A2(i1, i1) += k;
		//A2(i1, i2) += -k;
		//A2(i2, i1) += -k;
		//A2(i2, i2) += k;
		//A2 += k* A* A.transpose();
	}
	m_L.setFromTriplets(tripletList.begin(), tripletList.end());
	/*
	#pragma omp for schedule(static) 
	for (unsigned int i = 0; i < A2.rows(); i++)
		for (unsigned int j = 0; j < A2.cols(); j++)
		{
			m_L.block(i * I.rows(), j * I.cols(), I.rows(), I.cols()) = A2(i, j) * I;
		}
	*/
}

void Liu13_ClothModel::setJMatrix()
{
	Real k = 1;
	std::vector<Spring> edge = m_springs;
	std::vector<Triplet> tripletList;
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	m_J.resize(size * 3, edge.size() * 3);

	for (unsigned i = 0; i < edge.size(); i++)
	{
		unsigned int i1 = edge[i].m_vert[0];
		unsigned int i2 = edge[i].m_vert[1];
		k = edge[i].coeff;
		//A.setZero();
		//S.setZero();
		//A(i1) = 1;
		//A(i2) = -1;
		//S(i) = 1;
		tripletList.push_back(Triplet(3 * i1, 3 * i, k));
		tripletList.push_back(Triplet(3 * i1 + 1, 3 * i + 1, k));
		tripletList.push_back(Triplet(3 * i1 + 2, 3 * i + 2, k));

		tripletList.push_back(Triplet(3 * i2, 3 * i, -k));
		tripletList.push_back(Triplet(3 * i2 + 1, 3 * i + 1, -k));
		tripletList.push_back(Triplet(3 * i2 + 2, 3 * i + 2, -k));

		//AS(i1, i) += k;
		//AS(i2, i) += -k;
		//AS += k * A * S.transpose();
	}
	m_J.setFromTriplets(tripletList.begin(), tripletList.end());
	/*
	#pragma omp for schedule(static) 
	for (unsigned int i = 0; i < AS.rows(); i++)
		for (unsigned int j = 0; j < AS.cols(); j++)
		{
			m_J.block(i * I.rows(), j * I.cols(), I.rows(), I.cols()) = AS(i, j) * I;
		}
	*/

}

void Liu13_ClothModel::setMLMatrix(Real h)
{
	m_ML = m_massMatrix + h * h * m_L;
	solver.compute(m_ML);
	if (solver.info() != Eigen::Success)
		printf("fail\n");
}

void Liu13_ClothModel::setRestLength(ParticleData& pd)
{
	std::vector<Spring> edge = m_springs;

	for (unsigned int i = 0; i < edge.size(); i++)
	{
		const unsigned int v1 = edge[i].m_vert[0];
		const unsigned int v2 = edge[i].m_vert[1];
		m_restLength.push_back((pd.getPosition(v1) - pd.getPosition(v2)).norm());
	}
}

void Liu13_ClothModel::getForceVector(ParticleData& pd, Real h)
{
	long long m;
	m = size;
	Eigen::VectorXd y;
	y.resize(3 * m);
	m_b.setZero();
	for (long long i = 0; i < m; i++)
	{
		Vector3r fext = h * h * pd.getForce(i);
		Vector3r force = 1 * fext;
		pd.setForce(i, Vector3r(0,0,0));
		
		m_b(3 * i) = force.x();
		m_b(3 * i + 1) = force.y();
		m_b(3 * i + 2) = force.z();

	}
	double damp = 0.99;
	y = (1+damp) * m_X - damp * m_XOld;
	m_b += m_massMatrix * y;
	
	m_XOld = m_X;
}

void Liu13_ClothModel::setVectorSize()
{
	m_X.resize(3 * size);
	m_XOld.resize(3 * size);
	m_b.resize(3 * size);
	m_d.resize(3 * m_springs.size());
	m_ML.resize(3 * size, 3 * size);
	m_fixed_points.insert(0);
	m_fixed_points.insert(49);
}

void Liu13_ClothModel::localStep(ParticleData& pd)
{
	std::vector<Spring> edge = m_springs;
	m_d.setZero();
	
	for (unsigned int i = 0; i < edge.size(); i++)
	{
		long long v1 = edge[i].m_vert[0];
		long long v2 = edge[i].m_vert[1];
		Vector3r p12(m_X(3*v1) - m_X(3*v2), m_X(3*v1 + 1) - m_X(3*v2 + 1), m_X(3*v1 + 2) - m_X(3*v2 + 2));
		Real r = m_restLength[i];

		p12.normalize();
		Vector3r d = r * p12;
		m_d(3 * i) = d.x();
		m_d(3 * i + 1) = d.y();
		m_d(3 * i + 2) = d.z();
	}
}

void Liu13_ClothModel::globalStep(ParticleData& pd, Real h)
{
	long long m;
	m = size;
	Eigen::VectorXd b(m * 3);
	b = h * h * m_J * m_d + m_b;
	m_X = solver.solve(b);
	if (solver.info() != Eigen::Success)
		printf("fail2\n");

	//fixOverSpring(pd);
	//vectorToPosition(pd);
	/*
	#pragma omp parallel
	#pragma omp for schedule(static) 
	for (long long i = 0; i < m; i++)
	{

		if (i < 50) continue;
		Vector3r db(m_b(3 * i), m_b(3 * i + 1), m_b(3 * i + 2));
		Vector3r dJ(Jd(3 * i), Jd(3 * i + 1), Jd(3 * i + 2));
		Vector3r dml1 = m_X.transpose() * m_ML.middleCols(3 * i, 3);
		//Vector3r dml2 = m_ML.middleRows(3 * i, 3) * m_X;
		Matrix3r dm = m_ML.block(3 * i, 3 * i, 3, 3);
		Vector3r x = m_X.segment(3 * i, 3);
		Vector3r dmx1 = dm * x;
		//Vector3r dmx2 = x.transpose() * dm;
		Vector3r dg = dJ -db - dml1 + dmx1;
		pd.getLastPosition(i) = pd.getOldPosition(i);
		pd.getOldPosition(i) = pd.getPosition(i);
		pd.getPosition(i) = dm.inverse() * dg;
		
	}
	*/
	
}

void Liu13_ClothModel::vectorToPosition(ParticleData& pd)
{
	Eigen::VectorXd m_v = (m_X - m_XOld)/0.005;
	for (long long i = 0; i < size; i++)
	{
		pd.getLastPosition(i) = pd.getOldPosition(i);
		pd.getOldPosition(i) = pd.getPosition(i);
		pd.getPosition(i) = Vector3r(m_X(3 * i), m_X(3 * i + 1), m_X(3 * i + 2));
		//pd.getVelocity(i) = Vector3r(m_v(3 * i), m_v(3 * i + 1), m_v(3 * i + 2));
	}
}

void Liu13_ClothModel::setFixPoint(ParticleData& pd)
{
	for (unsigned int i = 0; i < size; i++)
	{
		if (m_fixed_points.find(i) != m_fixed_points.end())
		{
			Vector3r pos = pd.getPosition(i);
			m_X(3 * i) = pos.x();
			m_X(3 * i + 1) = pos.y();
			m_X(3 * i + 2) = pos.z();
		}
	}
}

void Liu13_ClothModel::fixOverSpring(ParticleData& pd, Real h2)
{
	std::vector<Spring> edge = m_springs;
	std::vector<unsigned int> save;
	for (unsigned int i = 0; i < edge.size(); i++)
	{
		if (edge[i].type != Springtype::st) continue;
		long long v1 = edge[i].m_vert[0];
		long long v2 = edge[i].m_vert[1];
		Vector3r p12(m_X(3 * v1) - m_X(3 * v2), m_X(3 * v1 + 1) - m_X(3 * v2 + 1), m_X(3 * v1 + 2) - m_X(3 * v2 + 2));
		Real r = m_restLength[i];
		Real diff = p12.norm() - r;
		if (diff > 0.1 * r)
		{
			p12.normalize();
			p12 *= 0.5 * (diff - 0.1 * r);
			/*
			if(v1 >= 50)
				m_X.segment(3 * v1, 3) -= 2* p12;
			if(v2 >= 50)
				m_X.segment(3 * v2, 3) += 2* p12;
			*/
			bool v1fix = (m_fixed_points.find(v1) != m_fixed_points.end());
			bool v2fix = (m_fixed_points.find(v2) != m_fixed_points.end());
			if (v1fix && v2fix) continue;
			else if(!v1fix && !v2fix)
			{
				pd.getForce(v1) -= p12 / h2;
				//m_X.segment(3 * v1, 3) = pd.getPosition(m_pindices[v1]);
				//save.push_back(v1);
				pd.getForce(v2) += p12 / h2;
				//m_X.segment(3 * v2, 3) = pd.getPosition(m_pindices[v2]);
				//save.push_back(v2);
			}
			else if (v2fix)
			{
				pd.getForce(v1) -= 2 * p12 / h2;
				//m_X.segment(3 * v1, 3) = pd.getPosition(m_pindices[v1]);
				//save.push_back(v1);
			}
			else 
			{
				pd.getForce(v2) += 2 * p12 / h2;
				//m_X.segment(3 * v2, 3) = pd.getPosition(m_pindices[v2]);
				//.push_back(v2);
			}

		}

	}
	/*
	for (unsigned int i = 0; i < save.size(); i++)
	{
		m_X.segment(3 * save[i], 3) = pd.getPosition(m_pindices[save[i]]);
	}
	*/
}

void Liu13_ClothModel::collisionSphere(ParticleData& pd)
{
	Vector3r center(4.5, -5, 4.5);
	Real radius = 3;
	for (unsigned int i = 0; i < size; i++)
	{
		Vector3r pos = pd.getPosition(i);
		if ((center - pos).norm() > radius) continue;
		//Vector3r oldp = m_XOld.segment(3 * i, 3);
		//Real length = (oldp - pos).norm();
		Vector3r dir = (pos - center).normalized();
		pd.getPosition(i) = center + dir * (radius+0.01);
		//m_XOld.segment(3 * i, 3) = pd.getPosition(i) - dir * length;
		//m_XOld.segment(3 * i, 3) = m_X.segment(3 * i, 3);
		m_X.segment(3 * i, 3) = pd.getPosition(i).cast<double>();
	}
}
void Liu13_ClothModel::restart()
{
	getPositionVector(*m_pd);
	m_XOld = m_X;
	printf("restart\n");
}
//
//void Liu13_ClothModel::fixedPointMovement(ParticleData& pd)
//{
//	static Real dx = -0.1;
//	bool b = false;
//	for (unsigned int i : m_fixed_points)
//	{
//		pd.getPosition(i).x() += dx;
//		if (abs(pd.getPosition(i).x()) > 3) b=true;
//	}
//	if (b) dx *= -1;
//}