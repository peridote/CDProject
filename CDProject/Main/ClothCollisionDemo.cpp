#include "Common/Common.h"
#include "Visualization/MiniGL.h"
#include "Visualization/Selection.h"
#include "GL/glut.h"
#include "Simulation/TimeManager.h"
#include <Eigen/Dense>
#include "Simulation/SimulationModel.h"
#include "Simulation/TimeStepController.h"
#include "Visualization/Visualization.h"
#include "Simulation/DistanceFieldCollisionDetection.h"
#include "Utils/OBJLoader.h"
#include "Utils/Logger.h"
#include "Utils/Timing.h"
#include "Utils/FileSystem.h"
#include "Common/DemoBase.h"
#include "Common/TweakBarParameters.h"
#include "Simulation/Simulation.h"
#include "TestManager/ImguiManager.h"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW 
#endif

using namespace PBD;
using namespace Eigen;
using namespace std;
using namespace Utilities;

void timeStep();
void buildModel();
void createMesh();
void render();
void reset();
void initParameters();
void exportMeshOBJ();
void exportOBJ();
void TW_CALL setBendingMethod(const void* value, void* clientData);
void TW_CALL getBendingMethod(void* value, void* clientData);
void TW_CALL setSimulationMethod(const void* value, void* clientData);
void TW_CALL getSimulationMethod(void* value, void* clientData);


const int nRows = 50;
const int nCols = 50;
const Real width = 10.0;
const Real height = 10.0;
short simulationMethod = 2;
short bendingMethod = 2;
bool doPause = true;
bool enableExportOBJ = false;
unsigned int exportFPS = 25;
Real nextFrameTime = 0.0;
unsigned int frameCounter = 1;
DemoBase* base;
DemoBase* base2;
DistanceFieldCollisionDetection cd;
//GLFWwindow* window = NULL; 

// main 
int main(int argc, char** argv)
{
	REPORT_MEMORY_LEAKS

	base = new DemoBase();
	base->init(argc, argv, "Cloth simulation");


	SimulationModel* model = new SimulationModel();
	model->init();
	Simulation::getCurrent()->setModel(model);

	buildModel();

	initParameters();

	Simulation::getCurrent()->setSimulationMethodChangedCallback([&]() { reset(); initParameters(); base->getSceneLoader()->readParameterObject(Simulation::getCurrent()->getTimeStep()); });

	//// OpenGL
	MiniGL::setClientIdleFunc(50, timeStep);
	MiniGL::setKeyFunc(0, 'r', reset);

	MiniGL::setClientSceneFunc(render);
	MiniGL::setViewport(40.0f, 0.1f, 500.0f, Vector3r(7.0, 4.0, 37.0), Vector3r(5.0, 0.0, 0.0));
	//MiniGL::setViewport(40.0f, 0.1f, 1000.0f, Vector3r(-10.0, 0.0, 5.0), Vector3r(0.0, 0.0, 0.0));

	TwType enumType2 = TwDefineEnum("SimulationMethodType", NULL, 0);
	TwAddVarCB(MiniGL::getTweakBar(), "SimulationMethod", enumType2, setSimulationMethod, getSimulationMethod, &simulationMethod, " label='Simulation method' enum='0 {None}, 1 {Distance constraints}, 2 {FEM based PBD}, 3 {Strain based dynamics}' group=Simulation");
	TwType enumType3 = TwDefineEnum("BendingMethodType", NULL, 0);
	TwAddVarCB(MiniGL::getTweakBar(), "BendingMethod", enumType3, setBendingMethod, getBendingMethod, &bendingMethod, " label='Bending method' enum='0 {None}, 1 {Dihedral angle}, 2 {Isometric bending}' group=Bending");

	ImguiManager* im = new ImguiManager();
	im->Initialize(base->getWindow());

	int my_image_width = 500;
	int my_image_height = 500;
	GLuint my_image_texture = 0;


	while (!glfwWindowShouldClose(im->m_window))
	{
		{

			im->StartFrame();
			im->createMainMenuBar();
			im->createLeftSideMenu();
			//im->createRightSideMenu();
			//im->createCenterMenu();
			im->createBottomMenu();
		}

		// render to texture
		im->fbo_bind();

		MiniGL::display();
		timeStep();

		im->fbo_unbind();


		// camera Menu
		static float trans_x = 7.0f;
		static float trans_y = 4.0f;
		static float trans_z = 37.0f;
		static float rot_w = 1.0f;
		static float rot_x = 0.0f;
		static float rot_y = 0.0f;
		static float rot_z = 0.0f;
		static float prev_rot_x = 0.0f;
		static float prev_rot_y = 0.0f;

		static int onOff = 0;
		ImGui::Begin("Controller");
		ImGui::SliderInt("on/off", &onOff, 0, 1);
		if (onOff == 0)
			base->m_doPause = true;
		else
			base->m_doPause = false;
		ImGui::Checkbox("exportObj", &enableExportOBJ);

		if (ImGui::Button("Reset"))
		{
			reset();
		}

		ImGui::SliderFloat("trans_x", &trans_x, -50.0f, 50.0f);
		ImGui::SliderFloat("trans_y", &trans_y, -50.0f, 50.0f);
		ImGui::SliderFloat("trans_z", &trans_z, -300.0f, 300.0f);
		MiniGL::m_translation.x() = -(Real)trans_x;
		MiniGL::m_translation.y() = -(Real)trans_y;
		MiniGL::m_translation.z() = -(Real)trans_z;

		//ImGui::SliderFloat("rot_w", &rot_w, 0.0f, 1.0f);
		ImGui::SliderFloat("rot_x", &rot_x, 0.0f, 1.0f);
		ImGui::SliderFloat("rot_y", &rot_y, 0.0f, 1.0f);
		if (rot_x != prev_rot_x) {
			if (rot_x < prev_rot_x)
				MiniGL::rotateX(-rot_x / static_cast<Real>(10.0));
			else
				MiniGL::rotateX(rot_x / static_cast<Real>(10.0));
			prev_rot_x = rot_x;
		}
		if (rot_y != prev_rot_y) {
			if (rot_y < prev_rot_y)
				MiniGL::rotateY(-rot_y / static_cast<Real>(10.0));
			else
				MiniGL::rotateY(rot_y / static_cast<Real>(10.0));
			prev_rot_y = rot_y;
		}

		ImGui::End();


		if (ImGui::Begin("Scene1"))
		{
			// dock layout by hard-coded or .ini file
			ImGui::BeginDockspace();

			if (ImGui::BeginDock("Scene_1")) {

				ImVec2 wsize = ImGui::GetWindowSize();
				MiniGL::width = im->m_w;
				MiniGL::height = im->m_h;
				ImGui::Image((void*)(intptr_t)im->m_fbo_texture, ImVec2(wsize.x, wsize.y), ImVec2(0, 1), ImVec2(1, 0));

			}
			ImGui::EndDock();

			if (ImGui::BeginDock("Scene_1_graph")) {
				float x_data[1000] = { 1, 2, 3, 4, 5 };
				float y_data[1000] = { 1, 1, 1, 1, 1 };


				//ImGui::BulletText("Move your mouse to change the data!");
				ImGui::BulletText("This assumes 60 FPS. Higher FPS requires larger buffer size.");
				static ImguiManager::ScrollingBuffer sdata1, sdata2;
				static ImguiManager::RollingBuffer   rdata1, rdata2;
				ImVec2 mouse = ImGui::GetMousePos();
				static float t = 0;
				if (onOff) {
					t += ImGui::GetIO().DeltaTime;

					sdata1.AddPoint(t, 1 / 100);
					rdata1.AddPoint(t, 1 / 100);
					sdata2.AddPoint(t, 1 / 100);
					rdata2.AddPoint(t, 1 / 100);

				}
				else {
					sdata1.AddPoint(t, 0);
					rdata1.AddPoint(t, 0);
					sdata2.AddPoint(t, 0);
					rdata2.AddPoint(t, 0);
				}
				static float history = 10.0f;
				ImGui::SliderFloat("History", &history, 1, 30, "%.1f s");
				rdata1.Span = history;
				rdata2.Span = history;

				static ImPlotAxisFlags flags = ImPlotAxisFlags_NoTickLabels;
				ImPlot::SetNextPlotLimitsX(t - history, t, ImGuiCond_Always);
				ImPlot::SetNextPlotLimitsY(0, 1);

				if (ImPlot::BeginPlot("##Scrolling", NULL, NULL, ImVec2(-1, 150), 0)) {
					ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
					ImPlot::PlotShaded("Relative Error", &sdata1.Data[0].x, &sdata1.Data[0].y, sdata1.Data.size(), -INFINITY, sdata1.Offset, 2 * sizeof(float));
					//ImPlot::PlotLine("Mouse Y", &sdata2.Data[0].x, &sdata2.Data[0].y, sdata2.Data.size(), sdata2.Offset, 2 * sizeof(float));
					ImPlot::EndPlot();
				}
				ImPlot::SetNextPlotLimitsX(0, history, ImGuiCond_Always);
				ImPlot::SetNextPlotLimitsY(0, 1);
				if (ImPlot::BeginPlot("##Rolling", NULL, NULL, ImVec2(-1, 150), 0, flags, flags)) {
					ImPlot::PlotLine("Relative Error", &rdata1.Data[0].x, &rdata1.Data[0].y, rdata1.Data.size(), 0, 2 * sizeof(float));
					//ImPlot::PlotLine("Mouse Y", &rdata2.Data[0].x, &rdata2.Data[0].y, rdata2.Data.size(), 0, 2 * sizeof(float));
					ImPlot::EndPlot();
				}


			}
			ImGui::EndDock();

			ImGui::EndDockspace();
		}

		if (ImGui::Begin("Scene2"))
		{
			// dock layout by hard-coded or .ini file
			ImGui::BeginDockspace();

			if (ImGui::BeginDock("Scene_2")) {

				ImVec2 wsize2 = ImGui::GetWindowSize();
				MiniGL::width = im->m_w;
				MiniGL::height = im->m_h;
				ImGui::Image((void*)(intptr_t)im->m_fbo_texture, ImVec2(wsize2.x, wsize2.y), ImVec2(0, 1), ImVec2(1, 0));

			}
			ImGui::EndDock();

			if (ImGui::BeginDock("Scene_2_graph")) {
				float x_data2[1000] = { 1, 2, 3, 4, 5 };
				float y_data2[1000] = { 1, 1, 1, 1, 1 };


				//ImGui::BulletText("Move your mouse to change the data!");
				ImGui::BulletText("This assumes 60 FPS. Higher FPS requires larger buffer size.");
				static ImguiManager::ScrollingBuffer sdata1, sdata2;
				static ImguiManager::RollingBuffer   rdata1, rdata2;
				ImVec2 mouse = ImGui::GetMousePos();
				static float t = 0;
				if (onOff) {
					t += ImGui::GetIO().DeltaTime;

					sdata1.AddPoint(t, 1 / 100);
					rdata1.AddPoint(t, 1 / 100);
					sdata2.AddPoint(t, 1 / 100);
					rdata2.AddPoint(t, 1 / 100);
				}
				else {
					sdata1.AddPoint(t, 0);
					rdata1.AddPoint(t, 0);
					sdata2.AddPoint(t, 0);
					rdata2.AddPoint(t, 0);
				}
				static float history = 10.0f;
				ImGui::SliderFloat("History", &history, 1, 30, "%.1f s");
				rdata1.Span = history;
				rdata2.Span = history;


				static ImPlotAxisFlags flags = ImPlotAxisFlags_NoTickLabels;
				ImPlot::SetNextPlotLimitsX(t - history, t, ImGuiCond_Always);
				ImPlot::SetNextPlotLimitsY(0, 1);

				if (ImPlot::BeginPlot("##Scrolling", NULL, NULL, ImVec2(-1, 150), 0)) {
					ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
					ImPlot::PlotShaded("Relative Error", &sdata1.Data[0].x, &sdata1.Data[0].y, sdata1.Data.size(), -INFINITY, sdata1.Offset, 2 * sizeof(float));
					//ImPlot::PlotLine("Mouse Y", &sdata2.Data[0].x, &sdata2.Data[0].y, sdata2.Data.size(), sdata2.Offset, 2 * sizeof(float));
					ImPlot::EndPlot();
				}
				ImPlot::SetNextPlotLimitsX(0, history, ImGuiCond_Always);
				ImPlot::SetNextPlotLimitsY(0, 1);
				if (ImPlot::BeginPlot("##Rolling", NULL, NULL, ImVec2(-1, 150), 0, flags, flags)) {
					ImPlot::PlotLine("Relative Error", &rdata1.Data[0].x, &rdata1.Data[0].y, rdata1.Data.size(), 0, 2 * sizeof(float));
					//ImPlot::PlotLine("Mouse Y", &rdata2.Data[0].x, &rdata2.Data[0].y, rdata2.Data.size(), 0, 2 * sizeof(float));
					ImPlot::EndPlot();
				}


			}
			ImGui::EndDock();

			ImGui::EndDockspace();
		}
		ImGui::End();

		/*static std::vector<float> v = { 1 };
		static float i = 0;
		v.push_back(i++);
		float* x_data = &v[0];
		static std::vector<float> w = { 3 };
		float* y_data = &w[0];*/

		im->Render();
	}
	im->Cleanup();
	delete(im);

	base->cleanup();

	Utilities::Timing::printAverageTimes();
	Utilities::Timing::printTimeSums();

	delete Simulation::getCurrent();
	delete base;
	delete model;

	return 0;
}


void initParameters()
{
	TwRemoveAllVars(MiniGL::getTweakBar());
	TweakBarParameters::cleanup();

	MiniGL::initTweakBarParameters();

	TweakBarParameters::createParameterGUI();
	TweakBarParameters::createParameterObjectGUI(base);
	TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent());
	TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getModel());
	TweakBarParameters::createParameterObjectGUI(Simulation::getCurrent()->getTimeStep());
}

void reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Simulation::getCurrent()->reset();
	base->getSelectedParticles().clear();

	Simulation::getCurrent()->getModel()->cleanup();
	Simulation::getCurrent()->getTimeStep()->getCollisionDetection()->cleanup();

	buildModel();
}

void timeStep()
{
	const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
	if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
		base->setValue(DemoBase::PAUSE, true);

	if (base->getValue<bool>(DemoBase::PAUSE))
		return;

	// Simulation code
	SimulationModel* model = Simulation::getCurrent()->getModel();
	const unsigned int numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
	for (unsigned int i = 0; i < numSteps; i++)
	{
		START_TIMING("SimStep");
		Simulation::getCurrent()->getTimeStep()->step(*model);
		STOP_TIMING_AVG;

		exportOBJ();
	}

	for (unsigned int i = 0; i < model->getTriangleModels().size(); i++)
		model->getTriangleModels()[i]->updateMeshNormals(model->getParticles());
}

void loadObj(const std::string& filename, VertexData& vd, IndexedFaceMesh& mesh, const Vector3r& scale)
{
	std::vector<OBJLoader::Vec3f> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<OBJLoader::Vec2f> texCoords;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	const unsigned int nTexCoords = (unsigned int)texCoords.size();
	mesh.initMesh(nPoints, nFaces * 2, nFaces);
	vd.reserve(nPoints);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		vd.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nTexCoords; i++)
	{
		mesh.addUV(texCoords[i][0], texCoords[i][1]);
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		// Reduce the indices by one
		int posIndices[3];
		int texIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j] - 1;
			if (nTexCoords > 0)
			{
				texIndices[j] = faces[i].texIndices[j] - 1;
				mesh.addUVIndex(texIndices[j]);
			}
		}

		mesh.addFace(&posIndices[0]);
	}
	mesh.buildNeighbors();

	mesh.updateNormals(vd, 0);
	mesh.updateVertexNormals(vd);

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
}

void buildModel()
{
	TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.005));

	createMesh();

	// create static rigid body
	string fileName = FileSystem::normalizePath(base->getDataPath() + "/models/cube.obj");
	IndexedFaceMesh mesh;
	VertexData vd;
	loadObj(fileName, vd, mesh, Vector3r::Ones());

	string fileNameTorus = FileSystem::normalizePath(base->getDataPath() + "/models/torus.obj");
	IndexedFaceMesh meshTorus;
	VertexData vdTorus;
	loadObj(fileNameTorus, vdTorus, meshTorus, Vector3r::Ones());

	SimulationModel* model = Simulation::getCurrent()->getModel();
	SimulationModel::RigidBodyVector& rb = model->getRigidBodies();

	rb.resize(2);
	// floor
	rb[0] = new RigidBody();
	rb[0]->initBody(1.0,
		Vector3r(0.0, -5.5, 0.0),
		Quaternionr(1.0, 0.0, 0.0, 0.0),
		vd, mesh,
		Vector3r(100.0, 1.0, 100.0));
	rb[0]->setMass(0.0);

	// torus
	rb[1] = new RigidBody();
	rb[1]->initBody(1.0,
		Vector3r(5.0, -1.5, 5.0),
		Quaternionr(1.0, 0.0, 0.0, 0.0),
		vdTorus, meshTorus,
		Vector3r(2.0, 2.0, 2.0));
	rb[1]->setMass(0.0);
	rb[1]->setFrictionCoeff(static_cast<Real>(0.1));

	Simulation::getCurrent()->getTimeStep()->setCollisionDetection(*model, &cd);
	cd.setTolerance(static_cast<Real>(0.05));

	const std::vector<Vector3r>* vertices1 = rb[0]->getGeometry().getVertexDataLocal().getVertices();
	const unsigned int nVert1 = static_cast<unsigned int>(vertices1->size());
	cd.addCollisionBox(0, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices1)[0], nVert1, Vector3r(100.0, 1.0, 100.0));

	const std::vector<Vector3r>* vertices2 = rb[1]->getGeometry().getVertexDataLocal().getVertices();
	const unsigned int nVert2 = static_cast<unsigned int>(vertices2->size());
	cd.addCollisionTorus(1, CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices2)[0], nVert2, Vector2r(2.0, 1.0));

	SimulationModel::TriangleModelVector& tm = model->getTriangleModels();
	ParticleData& pd = model->getParticles();
	for (unsigned int i = 0; i < tm.size(); i++)
	{
		const unsigned int nVert = tm[i]->getParticleMesh().numVertices();
		unsigned int offset = tm[i]->getIndexOffset();
		tm[i]->setFrictionCoeff(static_cast<Real>(0.1));
		cd.addCollisionObjectWithoutGeometry(i, CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
	}
}


void render()
{
	base->render();
}


/** Create a particle model mesh
*/
void createMesh()
{
	TriangleModel::ParticleMesh::UVs uvs;
	uvs.resize(nRows * nCols);

	const Real dy = width / (Real)(nCols - 1);
	const Real dx = height / (Real)(nRows - 1);

	Vector3r points[nRows * nCols];
	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			const Real y = (Real)dy * j;
			const Real x = (Real)dx * i;
			points[i * nCols + j] = Vector3r(x, 1.0, y);

			uvs[i * nCols + j][0] = x / width;
			uvs[i * nCols + j][1] = y / height;
		}
	}
	const int nIndices = 6 * (nRows - 1) * (nCols - 1);

	TriangleModel::ParticleMesh::UVIndices uvIndices;
	uvIndices.resize(nIndices);

	unsigned int indices[nIndices];
	int index = 0;
	for (int i = 0; i < nRows - 1; i++)
	{
		for (int j = 0; j < nCols - 1; j++)
		{
			int helper = 0;
			if (i % 2 == j % 2)
				helper = 1;

			indices[index] = i * nCols + j;
			indices[index + 1] = i * nCols + j + 1;
			indices[index + 2] = (i + 1) * nCols + j + helper;

			uvIndices[index] = i * nCols + j;
			uvIndices[index + 1] = i * nCols + j + 1;
			uvIndices[index + 2] = (i + 1) * nCols + j + helper;
			index += 3;

			indices[index] = (i + 1) * nCols + j + 1;
			indices[index + 1] = (i + 1) * nCols + j;
			indices[index + 2] = i * nCols + j + 1 - helper;

			uvIndices[index] = (i + 1) * nCols + j + 1;
			uvIndices[index + 1] = (i + 1) * nCols + j;
			uvIndices[index + 2] = i * nCols + j + 1 - helper;
			index += 3;
		}
	}

	SimulationModel* model = Simulation::getCurrent()->getModel();
	model->addTriangleModel(nRows * nCols, nIndices / 3, &points[0], &indices[0], uvIndices, uvs);

	ParticleData& pd = model->getParticles();
	for (unsigned int i = 0; i < pd.getNumberOfParticles(); i++)
	{
		pd.setMass(i, 1.0);
	}

	// init constraints
	for (unsigned int cm = 0; cm < model->getTriangleModels().size(); cm++)
	{
		if (simulationMethod == 1)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			const unsigned int nEdges = model->getTriangleModels()[cm]->getParticleMesh().numEdges();
			const IndexedFaceMesh::Edge* edges = model->getTriangleModels()[cm]->getParticleMesh().getEdges().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const unsigned int v1 = edges[i].m_vert[0] + offset;
				const unsigned int v2 = edges[i].m_vert[1] + offset;

				model->addDistanceConstraint(v1, v2);
			}
		}
		else if (simulationMethod == 2)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			TriangleModel::ParticleMesh& mesh = model->getTriangleModels()[cm]->getParticleMesh();
			const unsigned int* tris = mesh.getFaces().data();
			const unsigned int nFaces = mesh.numFaces();
			for (unsigned int i = 0; i < nFaces; i++)
			{
				const unsigned int v1 = tris[3 * i] + offset;
				const unsigned int v2 = tris[3 * i + 1] + offset;
				const unsigned int v3 = tris[3 * i + 2] + offset;
				model->addFEMTriangleConstraint(v1, v2, v3);
			}
		}
		else if (simulationMethod == 3)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			TriangleModel::ParticleMesh& mesh = model->getTriangleModels()[cm]->getParticleMesh();
			const unsigned int* tris = mesh.getFaces().data();
			const unsigned int nFaces = mesh.numFaces();
			for (unsigned int i = 0; i < nFaces; i++)
			{
				const unsigned int v1 = tris[3 * i] + offset;
				const unsigned int v2 = tris[3 * i + 1] + offset;
				const unsigned int v3 = tris[3 * i + 2] + offset;
				model->addStrainTriangleConstraint(v1, v2, v3);
			}
		}
		if (bendingMethod != 0)
		{
			const unsigned int offset = model->getTriangleModels()[cm]->getIndexOffset();
			TriangleModel::ParticleMesh& mesh = model->getTriangleModels()[cm]->getParticleMesh();
			unsigned int nEdges = mesh.numEdges();
			const TriangleModel::ParticleMesh::Edge* edges = mesh.getEdges().data();
			const unsigned int* tris = mesh.getFaces().data();
			for (unsigned int i = 0; i < nEdges; i++)
			{
				const int tri1 = edges[i].m_face[0];
				const int tri2 = edges[i].m_face[1];
				if ((tri1 != 0xffffffff) && (tri2 != 0xffffffff))
				{
					// Find the triangle points which do not lie on the axis
					const int axisPoint1 = edges[i].m_vert[0];
					const int axisPoint2 = edges[i].m_vert[1];
					int point1 = -1;
					int point2 = -1;
					for (int j = 0; j < 3; j++)
					{
						if ((tris[3 * tri1 + j] != axisPoint1) && (tris[3 * tri1 + j] != axisPoint2))
						{
							point1 = tris[3 * tri1 + j];
							break;
						}
					}
					for (int j = 0; j < 3; j++)
					{
						if ((tris[3 * tri2 + j] != axisPoint1) && (tris[3 * tri2 + j] != axisPoint2))
						{
							point2 = tris[3 * tri2 + j];
							break;
						}
					}
					if ((point1 != -1) && (point2 != -1))
					{
						const unsigned int vertex1 = point1 + offset;
						const unsigned int vertex2 = point2 + offset;
						const unsigned int vertex3 = edges[i].m_vert[0] + offset;
						const unsigned int vertex4 = edges[i].m_vert[1] + offset;
						if (bendingMethod == 1)
							model->addDihedralConstraint(vertex1, vertex2, vertex3, vertex4);
						else if (bendingMethod == 2)
							model->addIsometricBendingConstraint(vertex1, vertex2, vertex3, vertex4);
					}
				}
			}
		}
	}

	LOG_INFO << "Number of triangles: " << nIndices / 3;
	LOG_INFO << "Number of vertices: " << nRows * nCols;

}

void exportMeshOBJ(const std::string& exportFileName, const unsigned int nVert, const Vector3r* pos, const unsigned int nTri, const unsigned int* faces)
{
	// Open the file
	std::ofstream outfile(exportFileName);
	if (!outfile)
	{
		LOG_WARN << "Cannot open a file to save OBJ mesh.";
		return;
	}

	// Header
	outfile << "# Created by the PositionBasedDynamics library\n";
	outfile << "g default\n";

	// Vertices
	{
		for (auto j = 0u; j < nVert; j++)
		{
			const Vector3r& x = pos[j];
			outfile << "v " << x[0] << " " << x[1] << " " << x[2] << "\n";
		}
	}

	// faces
	{
		for (auto j = 0; j < nTri; j++)
		{
			outfile << "f " << faces[3 * j + 0] + 1 << " " << faces[3 * j + 1] + 1 << " " << faces[3 * j + 2] + 1 << "\n";
		}
	}
	outfile.close();
}

void exportOBJ()
{
	if (!enableExportOBJ)
		return;

	if (TimeManager::getCurrent()->getTime() < nextFrameTime)
		return;

	nextFrameTime += 1.0 / (Real)exportFPS;

	//////////////////////////////////////////////////////////////////////////
	// rigid bodies
	//////////////////////////////////////////////////////////////////////////

	std::string exportPath = base->getOutputPath() + "/export";
	FileSystem::makeDirs(exportPath);

	SimulationModel* model = Simulation::getCurrent()->getModel();
	const ParticleData& pd = model->getParticles();
	for (unsigned int i = 0; i < model->getTriangleModels().size(); i++)
	{
		const IndexedFaceMesh& mesh = model->getTriangleModels()[i]->getParticleMesh();
		const unsigned int offset = model->getTriangleModels()[i]->getIndexOffset();
		const Vector3r* x = model->getParticles().getVertices()->data();

		std::string fileName = "triangle_model";
		fileName = fileName + std::to_string(i) + "_" + std::to_string(frameCounter) + ".obj";
		std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

		exportMeshOBJ(exportFileName, mesh.numVertices(), &x[offset], mesh.numFaces(), mesh.getFaces().data());
	}

	for (unsigned int i = 0; i < model->getTetModels().size(); i++)
	{
		const IndexedFaceMesh& mesh = model->getTetModels()[i]->getVisMesh();
		const unsigned int offset = model->getTetModels()[i]->getIndexOffset();
		const Vector3r* x = model->getTetModels()[i]->getVisVertices().getVertices()->data();

		std::string fileName = "tet_model";
		fileName = fileName + std::to_string(i) + "_" + std::to_string(frameCounter) + ".obj";
		std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

		exportMeshOBJ(exportFileName, mesh.numVertices(), x, mesh.numFaces(), mesh.getFaces().data());
	}

	for (unsigned int i = 0; i < model->getRigidBodies().size(); i++)
	{
		const IndexedFaceMesh& mesh = model->getRigidBodies()[i]->getGeometry().getMesh();
		const Vector3r* x = model->getRigidBodies()[i]->getGeometry().getVertexData().getVertices()->data();

		std::string fileName = "rigid_body";
		fileName = fileName + std::to_string(i) + "_" + std::to_string(frameCounter) + ".obj";
		std::string exportFileName = FileSystem::normalizePath(exportPath + "/" + fileName);

		exportMeshOBJ(exportFileName, mesh.numVertices(), x, mesh.numFaces(), mesh.getFaces().data());
	}


	frameCounter++;
}

void TW_CALL setBendingMethod(const void* value, void* clientData)
{
	const short val = *(const short*)(value);
	*((short*)clientData) = val;
	reset();
}

void TW_CALL getBendingMethod(void* value, void* clientData)
{
	*(short*)(value) = *((short*)clientData);
}

void TW_CALL setSimulationMethod(const void* value, void* clientData)
{
	const short val = *(const short*)(value);
	*((short*)clientData) = val;
	reset();
}

void TW_CALL getSimulationMethod(void* value, void* clientData)
{
	*(short*)(value) = *((short*)clientData);
}