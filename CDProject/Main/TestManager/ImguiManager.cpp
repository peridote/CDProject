#include "Utils/OBJLoader.h"
#include "Simulation/SimulationModel.h"
#include "Utils/FileSystem.h"
#include "Simulation/Simulation.h"
#include "math.h"


#include "ImguiManager.h"

const int nRows = 50;
const int nCols = 50;
const Real width = 10.0;
const Real height = 10.0;
short simulationMethod = 2;
short bendingMethod = 2;
int nCloth = 0;
int nTorus = 0;
int nFloor = 0;
int nSphere = 0;
int nCylinder = 0;
int nCube = 0;
int nCustom = 0;
int direction = 0;

ImguiManager::ImguiManager()
{
	m_window = NULL;
	m_w = m_h = 0;
	m_clear_color = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
	m_filetreeitem_current_idx = 0;
	m_filetree_double_clicked_item = "";
	//m_torus_tree_items = {""};
	m_rigidbody_num = 1; // ground가 이미 있어서 1로 초기화
	which_tree = 6;

	translation = Eigen::Vector3f(0, 0, 0);

	window_flags = 0;
	//if (no_titlebar)        window_flags |= ImGuiWindowFlags_NoTitleBar;
	//if (no_scrollbar)       window_flags |= ImGuiWindowFlags_NoScrollbar;
	//if (!no_menu)           window_flags |= ImGuiWindowFlags_MenuBar;
	window_flags |= ImGuiWindowFlags_NoMove;
	window_flags |= ImGuiWindowFlags_NoResize;
	window_flags |= ImGuiWindowFlags_NoCollapse;
	//if (no_nav)             window_flags |= ImGuiWindowFlags_NoNav;
	//if (no_background)      window_flags |= ImGuiWindowFlags_NoBackground;
	//if (no_bring_to_front)  window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus;
	//if (no_close)           p_open = NULL; // Don't pass our bool* to Begin
}

ImguiManager::~ImguiManager()
{

}


void ImguiManager::draw()
{

}

void ImguiManager::Initialize(GLFWwindow* window)
{
	m_window = window;
	// Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
	// GL ES 2.0 + GLSL 100
	const char* glsl_version = "#version 100";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
	// GL 3.2 + GLSL 150
	const char* glsl_version = "#version 150";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
	// GL 3.0 + GLSL 130
	const char* glsl_version = "#version 130";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
#endif


	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImPlot::CreateContext();
	m_io = &ImGui::GetIO(); (void)m_io;
	//io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
	//io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsClassic();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(m_window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	// Load Fonts
	m_io->Fonts->AddFontDefault();
	m_io->Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
	m_io->Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
	m_io->Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
	m_io->Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
	ImFont* font = m_io->Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, m_io->Fonts->GetGlyphRangesJapanese());

	// Our state
	//m_clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
	m_clear_color = ImVec4(0.0f, 0.0f, 0.0f, 1.00f);

	m_fileDialog.SetTitle("title");
	m_fileDialog.SetTypeFilters({ ".h", ".cpp", ".obj" });

	glfwGetWindowSize(m_window, &m_w, &m_h);

	// Initialize OpenGL loader
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
	bool err = gl3wInit() != 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
	bool err = glewInit() != GLEW_OK;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
	bool err = gladLoadGL() == 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD2)
	bool err = gladLoadGL(glfwGetProcAddress) == 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING2)
	bool err = false;
	glbinding::Binding::initialize();
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING3)
	bool err = false;
	glbinding::initialize([](const char* name) { return (glbinding::ProcAddress)glfwGetProcAddress(name); });
#else
	bool err = false;
#endif
	if (err)
	{
		fprintf(stderr, "Failed to initialize OpenGL loader!\n");
		return exit(1);
	}

}

void ImguiManager::StartFrame()
{
	glfwPollEvents();
	//glfwSetKeyCallback(m_window, keyCallback);

	// Start the Dear ImGui frame
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();
}

void ImguiManager::Render()
{
	ImGui::Render();
	int display_w, display_h;
	glfwGetFramebufferSize(m_window, &display_w, &display_h);
	glViewport(0, 0, display_w, display_h);
	glClearColor(m_clear_color.x * m_clear_color.w, m_clear_color.y * m_clear_color.w, m_clear_color.z * m_clear_color.w, m_clear_color.w);
	glClear(GL_COLOR_BUFFER_BIT);
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	glfwSwapBuffers(m_window);
}

void ImguiManager::Cleanup()
{
	fbo_cleanup();
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImPlot::DestroyContext();
	ImGui::DestroyContext();

	glfwDestroyWindow(m_window);
	glfwTerminate();
}


void ImguiManager::reset()
{
	m_cube_tree_items.clear();
	m_cylinder_tree_items.clear();
	m_sphere_tree_items.clear();
	m_cloth_tree_items.clear();
	m_torus_tree_items.clear();
	m_floor_tree_items.clear();
	m_custom_tree_items.clear();
	m_map_filetree_rigidbody.clear();
	//m_custom_tree_items.push_back("(1)torus.obj"); // 현재 collision demo 기본 obj가 tours라서 넣어줌...
	m_rigidbody_num = 1;

	nCloth = 0;
	nTorus = 0;
	nFloor = 0;
	nSphere = 0;
	nCylinder = 0;
	nCube = 0;
	nCustom = 0;
}

void ImguiManager::createMainMenuBar()
{
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("New"))
			{
				//Do something
			}
			ImGui::EndMenu();
		}

		ImGui::EndMainMenuBar();
	}
}

/** Create a particle model mesh
*/
void ImguiManager::createMesh()
{
	PBD::TriangleModel::ParticleMesh::UVs uvs;
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

	PBD::TriangleModel::ParticleMesh::UVIndices uvIndices;
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

	PBD::SimulationModel* model = PBD::Simulation::getCurrent()->getModel();
	model->addTriangleModel(nRows * nCols, nIndices / 3, &points[0], &indices[0], uvIndices, uvs);

	PBD::ParticleData& pd = model->getParticles();
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
			const Utilities::IndexedFaceMesh::Edge* edges = model->getTriangleModels()[cm]->getParticleMesh().getEdges().data();
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
			PBD::TriangleModel::ParticleMesh& mesh = model->getTriangleModels()[cm]->getParticleMesh();
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
			PBD::TriangleModel::ParticleMesh& mesh = model->getTriangleModels()[cm]->getParticleMesh();
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
			PBD::TriangleModel::ParticleMesh& mesh = model->getTriangleModels()[cm]->getParticleMesh();
			unsigned int nEdges = mesh.numEdges();
			const PBD::TriangleModel::ParticleMesh::Edge* edges = mesh.getEdges().data();
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

	PBD::SimulationModel::TriangleModelVector& tm = model->getTriangleModels();
	
	for (unsigned int i = 0; i < tm.size(); i++)
	{
		const unsigned int nVert = tm[i]->getParticleMesh().numVertices();
		unsigned int offset = tm[i]->getIndexOffset();
		tm[i]->setFrictionCoeff(static_cast<Real>(0.1));
		PBD::Simulation* current = PBD::Simulation::getCurrent();
		if (current == PBD::Simulation::m_simul1) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep()->getCollisionDetection();
			cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
		} else if (current == PBD::Simulation::m_simul2) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep2()->getCollisionDetection();
			cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
		}
		//cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
	}

	LOG_INFO << "Number of triangles: " << nIndices / 3;
	LOG_INFO << "Number of vertices: " << nRows * nCols;

}

void ImguiManager::createLiuMesh()
{
	PBD::TriangleModel::ParticleMesh::UVs uvs;
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

	PBD::TriangleModel::ParticleMesh::UVIndices uvIndices;
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

	PBD::SimulationModel* model = PBD::Simulation::getCurrent()->getModel();
	model->addLiuModel(nRows * nCols, nIndices / 3, &points[0], &indices[0], uvIndices, uvs);
	ParticleData& pd = model->getParticles();
	Liu13_ClothModel* tri;
	
	tri = (Liu13_ClothModel*)(model->getTriangleModels()[0]);
	tri->init(pd);
	//
	PBD::SimulationModel::TriangleModelVector& tm = model->getTriangleModels();
	PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)PBD::Simulation::getCurrent()->getTimeStep()->getCollisionDetection();
	for (unsigned int i = 0; i < tm.size(); i++)
	{
		const unsigned int nVert = tm[i]->getParticleMesh().numVertices();
		unsigned int offset = tm[i]->getIndexOffset();
		tm[i]->setFrictionCoeff(static_cast<Real>(0.1));
		PBD::Simulation* current = PBD::Simulation::getCurrent();
		if (current == PBD::Simulation::m_simul1) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep()->getCollisionDetection();
			cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
		}
		else if (current == PBD::Simulation::m_simul2) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep2()->getCollisionDetection();
			cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
		}
		//cd.addCollisionObjectWithoutGeometry(i, PBD::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType, &pd.getPosition(offset), nVert, true);
	}

	//LOG_INFO << "Number of triangles: " << nIndices / 3;
	//LOG_INFO << "Number of vertices: " << nRows * nCols;

}

void ImguiManager::createTorus()
{
	addRigidbody(torus);
}

void ImguiManager::createLeftSideMenu()
{
	if (ImGui::Begin("Scene Hierarchy", NULL, window_flags))
	{
		ImGui::SetWindowPos(ImVec2(0, 20));
		ImGui::SetWindowSize(ImVec2(180, 837), 0);

		static ImGuiTreeNodeFlags base_flags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick | ImGuiTreeNodeFlags_SpanAvailWidth;
		static bool align_label_with_current_x_position = false;
		static bool test_drag_and_drop = false;

		ImGui::SetNextItemOpen(true, ImGuiCond_Once);
		if (ImGui::TreeNode("Cube"))
		{
			ImGui::SameLine(120, 0);
			if (ImGui::SmallButton("create"))
			{
				if (nCube < 10)
				{
					addRigidbody(cube);
					PBD::Simulation::switchCurrent();
					addRigidbody(cube);
					PBD::Simulation::switchCurrent();

					nCube += 1;
					std::string name = "(" + std::to_string(nCube) + ")" + "cube";
					
					m_map_filetree_rigidbody.insert(std::make_pair(name, m_rigidbody_num++));
					m_cube_tree_items.push_back(name);
				}
				else
				{
					// message 출력
				}
			}

			for (int n = 0; n < m_cube_tree_items.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n && which_tree == cube);
				if (ImGui::Selectable(m_cube_tree_items[n].c_str(), is_selected))
				{
					m_filetree_current_item = m_cube_tree_items[n];
					m_filetreeitem_current_idx = n;
					which_tree = cube;
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
			ImGui::TreePop();
		}

		ImGui::SetNextItemOpen(true, ImGuiCond_Once);
		if (ImGui::TreeNode("Cylinder"))
		{
			ImGui::SameLine(120,0);
			if (ImGui::SmallButton("create"))
			{
				if (nCylinder < 10)
				{
					addRigidbody(cylinder);
					PBD::Simulation::switchCurrent();
					addRigidbody(cylinder);
					PBD::Simulation::switchCurrent();

					nCylinder += 1;
					std::string name = "(" + std::to_string(nCylinder) + ")" + "cylinder";

					m_map_filetree_rigidbody.insert(std::make_pair(name, m_rigidbody_num++));
					m_cylinder_tree_items.push_back(name);
				}
				else
				{
					// message 출력
				}
			}

			for (int n = 0; n < m_cylinder_tree_items.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n && which_tree == cylinder);
				if (ImGui::Selectable(m_cylinder_tree_items[n].c_str(), is_selected))
				{
					m_filetree_current_item = m_cylinder_tree_items[n];
					m_filetreeitem_current_idx = n;
					which_tree = cylinder;
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
			ImGui::TreePop();
		}

		ImGui::SetNextItemOpen(true, ImGuiCond_Once);
		if (ImGui::TreeNode("Sphere"))
		{
			ImGui::SameLine(120, 0);
			if (ImGui::SmallButton("create")) 
			{
				if (nSphere < 10)
				{
					addRigidbody(sphere);
					PBD::Simulation::switchCurrent();
					addRigidbody(sphere);
					PBD::Simulation::switchCurrent();

					nSphere += 1;
					std::string name = "(" + std::to_string(nSphere) + ")" + "sphere";

					m_map_filetree_rigidbody.insert(std::make_pair(name, m_rigidbody_num++));
					m_sphere_tree_items.push_back(name);
				}
				else
				{
					// message 출력
				}
			}

			for (int n = 0; n < m_sphere_tree_items.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n && which_tree == sphere);
				if (ImGui::Selectable(m_sphere_tree_items[n].c_str(), is_selected))
				{
					m_filetree_current_item = m_sphere_tree_items[n];
					m_filetreeitem_current_idx = n;
					which_tree = sphere;
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
			ImGui::TreePop();
		}

		ImGui::SetNextItemOpen(true, ImGuiCond_Once);
		if (ImGui::TreeNode("Cloth"))
		{
			ImGui::SameLine(120, 0);
			if (ImGui::SmallButton("create")) 
			{
				if (nCloth == 0)
				{
					createMesh();
					PBD::Simulation::switchCurrent();
					SimulationMethods method = static_cast<SimulationMethods>(PBD::Simulation::getCurrent()->getSimulationMethod());

					if (method == SimulationMethods::PBD)
					{
						createMesh();
					}
					else if (method == SimulationMethods::TEST)
					{
						createLiuMesh();
					}
					
					PBD::Simulation::switchCurrent();

					nCloth += 1;
					m_cloth_tree_items.push_back("(" + std::to_string(nCloth) + ")" + "cloth");
				}
				else
				{
					// message 출력
				}
			}

			for (int n = 0; n < m_cloth_tree_items.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n && which_tree == cloth);
				if (ImGui::Selectable(m_cloth_tree_items[n].c_str(), is_selected))
				{
					m_filetree_current_item = m_cloth_tree_items[n];
					m_filetreeitem_current_idx = n;
					which_tree = cloth;
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
			ImGui::TreePop();
		}

		ImGui::SetNextItemOpen(true, ImGuiCond_Once);
		if (ImGui::TreeNode("Torus"))
		{
			ImGui::SameLine(120, 0);
			if (ImGui::SmallButton("create"))
			{
				if (nTorus == 0)
				{
				
					addRigidbody(torus);
					PBD::Simulation::switchCurrent();
					addRigidbody(torus);
					PBD::Simulation::switchCurrent();

					nTorus += 1;
					std::string name = "(" + std::to_string(nTorus) + ")" + "torus";

					m_map_filetree_rigidbody.insert(std::make_pair(name, m_rigidbody_num++));
					m_torus_tree_items.push_back(name);
				}
				else
				{
					// message 출력

				}
			}

			for (int n = 0; n < m_torus_tree_items.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n && which_tree == torus);
				if (ImGui::Selectable(m_torus_tree_items[n].c_str(), is_selected))
				{
					m_filetree_current_item = m_torus_tree_items[n];
					m_filetreeitem_current_idx = n;
					which_tree = torus;
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
			ImGui::TreePop();
		}

		ImGui::SetNextItemOpen(true, ImGuiCond_Once);
		if (ImGui::TreeNode("Floor"))
		{
			ImGui::SameLine(120, 0);
			if (ImGui::SmallButton("create"))
			{
				if (nFloor == 0)
				{

					addRigidbody(floor);
					PBD::Simulation::switchCurrent();
					addRigidbody(floor);
					PBD::Simulation::switchCurrent();

					nFloor += 1;
					std::string name = "(" + std::to_string(nTorus) + ")" + "floor";

					m_map_filetree_rigidbody.insert(std::make_pair(name, m_rigidbody_num++));
					m_floor_tree_items.push_back(name);
				}
				else
				{
					// message 출력

				}
			}

			for (int n = 0; n < m_floor_tree_items.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n && which_tree == floor);
				if (ImGui::Selectable(m_floor_tree_items[n].c_str(), is_selected))
				{
					m_filetree_current_item = m_floor_tree_items[n];
					m_filetreeitem_current_idx = n;
					which_tree = floor;
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
			ImGui::TreePop();
		}

		ImGui::SetNextItemOpen(true, ImGuiCond_Once);
		if (ImGui::TreeNode("Custom"))
		{
			if (m_filetree_double_clicked_item != "")
			{
				nCustom += 1;
				m_filetree_double_clicked_item = "(" + std::to_string(nCustom) + ")" + m_filetree_double_clicked_item;
				m_map_filetree_rigidbody.insert(std::make_pair(m_filetree_double_clicked_item, m_rigidbody_num++));
				m_custom_tree_items.push_back(m_filetree_double_clicked_item);
				std::cout << m_filetree_double_clicked_item << std::endl;
				m_filetree_double_clicked_item = "";
			}

			for (int n = 0; n < m_custom_tree_items.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n && which_tree == custom);
				if (ImGui::Selectable(m_custom_tree_items[n].c_str(), is_selected)) 
				{
					m_filetree_current_item = m_custom_tree_items[n];
					m_filetreeitem_current_idx = n;
					which_tree = custom;
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
			ImGui::TreePop();
		}
		//if (ImGui::TreeNode("Cube"))
		//{

		//	if (align_label_with_current_x_position)
		//		ImGui::Unindent(ImGui::GetTreeNodeToLabelSpacing());

		//	// 'selection_mask' is dumb representation of what may be user-side selection state.
		//	//  You may retain selection state inside or outside your objects in whatever format you see fit.
		//	// 'node_clicked' is temporary storage of what node we have clicked to process selection at the end
		//	/// of the loop. May be a pointer to your own node type, etc.
		//	static int selection_mask = (1 << 2);
		//	int node_clicked = -1;
		//	for (int i = 0; i < 6; i++)
		//	{
		//		// Disable the default "open on single-click behavior" + set Selected flag according to our selection.
		//		ImGuiTreeNodeFlags node_flags = base_flags;
		//		const bool is_selected = (selection_mask & (1 << i)) != 0;
		//		if (is_selected)
		//			node_flags |= ImGuiTreeNodeFlags_Selected;
		//		if (i < 3)
		//		{
		//			// Items 0..2 are Tree Node
		//			bool node_open = ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "Selectable Node %d", i);
		//			if (ImGui::IsItemClicked())
		//				node_clicked = i;
		//			if (test_drag_and_drop && ImGui::BeginDragDropSource())
		//			{
		//				ImGui::SetDragDropPayload("_TREENODE", NULL, 0);
		//				ImGui::Text("This is a drag and drop source");
		//				ImGui::EndDragDropSource();
		//			}
		//			if (node_open)
		//			{
		//				ImGui::BulletText("Blah blah\nBlah Blah");
		//				ImGui::TreePop();
		//			}
		//		}
		//		else
		//		{
		//			// Items 3..5 are Tree Leaves
		//			// The only reason we use TreeNode at all is to allow selection of the leaf. Otherwise we can
		//			// use BulletText() or advance the cursor by GetTreeNodeToLabelSpacing() and call Text().
		//			node_flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen; // ImGuiTreeNodeFlags_Bullet
		//			ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "Selectable Leaf %d", i);
		//			if (ImGui::IsItemClicked())
		//				node_clicked = i;
		//			if (test_drag_and_drop && ImGui::BeginDragDropSource())
		//			{
		//				ImGui::SetDragDropPayload("_TREENODE", NULL, 0);
		//				ImGui::Text("This is a drag and drop source");
		//				ImGui::EndDragDropSource();
		//			}
		//		}
		//	}
		//	if (node_clicked != -1)
		//	{
		//		// Update selection state
		//		// (process outside of tree loop to avoid visual inconsistencies during the clicking frame)
		//		if (ImGui::GetIO().KeyCtrl)
		//			selection_mask ^= (1 << node_clicked);          // CTRL+click to toggle
		//		else //if (!(selection_mask & (1 << node_clicked))) // Depending on selection behavior you want, may want to preserve selection when clicking on item that is part of the selection
		//			selection_mask = (1 << node_clicked);           // Click to single-select
		//	}
		//	if (align_label_with_current_x_position)
		//		ImGui::Indent(ImGui::GetTreeNodeToLabelSpacing());
		//	ImGui::TreePop();
		//}
		//ImGui::TreePop();
		// left Menu
		//if (ImGui::BeginDock("Project Files"))
	
		//ImGui::SetWindowPos(ImVec2(0, 20));
		//ImGui::SetWindowSize(ImVec2(m_w / 6, m_h / 3 * 2));
	}
	ImGui::End();
}

void ImguiManager::createRightSideMenu()
{

	if (ImGui::Begin("Properties", NULL, window_flags))
	{
		ImGui::SetWindowPos(ImVec2(180 + 2 * m_w / 2.5, 20 + 467));
		ImGui::SetWindowSize(ImVec2(204, 370));
		//ImGui::SetWindowSize(ImVec2(m_w / 6, m_h / 3 * 2));
		//ImGui::SetWindowPos(ImVec2(m_w / 6 * 4, 20));
		//ImGui::SetWindowSize(ImVec2(m_w / 6 * 2, m_h / 3 * 2), 0);
		if (ImGui::BeginTabBar("right"))
		{

			ImGui::PushItemWidth(ImGui::GetFontSize() * 8);
			if (strcmp(m_filetree_current_item.c_str(), "(1)torus") == 0)
			{
				static bool test1 = true;
				ImGui::Checkbox("Collision", &test1);
				static bool test2 = true;
				ImGui::Checkbox("Visible", &test2);
				static bool test3 = true;
				ImGui::Checkbox("Static", &test3);
				ImGui::EndTabBar();
			}
			else if (strcmp(m_filetree_current_item.c_str(), "(1)cloth") == 0)
			{
				static bool test1 = true;
				ImGui::Checkbox("Collision", &test1);
				static bool test2 = true;
				ImGui::Checkbox("Visible", &test2);
				static bool test3 = true;
				ImGui::Checkbox("Static", &test3);
				ImGui::EndTabBar();
			}
			else if (std::filesystem::path(m_filetree_current_item).extension() == ".obj" || (which_tree >= 0 && which_tree <= 2))
			{
				int n = m_map_filetree_rigidbody.find(m_filetree_current_item)->second;
				static bool test1 = true;
				ImGui::Checkbox("Collision", &test1);
				static bool test2 = true;
				ImGui::Checkbox("Visible", &test2);
				static bool test3 = true;

				PBD::SimulationModel* model = PBD::Simulation::getCurrent()->getModel();
				PBD::SimulationModel::RigidBodyVector& rb = model->getRigidBodies();

				static float rb_scale_x = 1.0f;
				static float rb_scale_y = 1.0f;
				static float rb_scale_z = 1.0f;
				static float prev_rb_scale_x = 1.0f;
				static float prev_rb_scale_y = 1.0f;
				static float prev_rb_scale_z = 1.0f;
				
				static float rb_trans_x = rb[n-1]->getPosition()(0);
				static float rb_trans_y = rb[n-1]->getPosition()(1);
				static float rb_trans_z = rb[n-1]->getPosition()(2);
				static float prev_rb_trans_x = 0.0f;
				static float prev_rb_trans_y = 0.0f;
				static float prev_rb_trans_z = 0.0f;

				static float rb_rot_w = rb[n-1]->getRotation().w();
				static float rb_rot_x = rb[n-1]->getRotation().x();
				static float rb_rot_y = rb[n-1]->getRotation().y();
				static float rb_rot_z = rb[n-1]->getRotation().z();
				static float prev_rb_rot_w = 1.0f;
				static float prev_rb_rot_x = 0.0f;
				static float prev_rb_rot_y = 0.0f;
				static float prev_rb_rot_z = 0.0f;

				for (int l = 0; l < IM_ARRAYSIZE(ImGui::GetIO().KeysDown); l++) {
					if (l == ImGui::GetKeyIndex(ImGuiKey_LeftArrow) && ImGui::IsKeyPressed(l))
					{
						printf("left\n");
						rb_trans_x -= 0.5f;
					}
					else if (l == ImGui::GetKeyIndex(ImGuiKey_RightArrow) && ImGui::IsKeyPressed(l))
					{
						printf("right\n");
						rb_trans_x += 0.5f;

					}
					else if (l == ImGui::GetKeyIndex(ImGuiKey_UpArrow) && ImGui::IsKeyPressed(l))
					{
						printf("up\n");
						rb_trans_y += 0.5f;
					}
					else if (l == ImGui::GetKeyIndex(ImGuiKey_DownArrow) && ImGui::IsKeyPressed(l))
					{
						printf("down\n");
						rb_trans_y -= 0.5f;
					}
				}

				ImGui::SliderFloat("rb_trans_x", &rb_trans_x, -5.0f, 5.0f);
				ImGui::SliderFloat("rb_trans_y", &rb_trans_y, -5.0f, 5.0f);
				ImGui::SliderFloat("rb_trans_z", &rb_trans_z, -5.0f, 5.0f);

				//ImGui::SliderFloat("rb_rot_w", &rb_rot_w, -1.0f, 1.0f);
				ImGui::SliderFloat("rb_rot_x", &rb_rot_x, 0.0f, 360.0f);
				ImGui::SliderFloat("rb_rot_y", &rb_rot_y, 0.0f, 360.0f);
				ImGui::SliderFloat("rb_rot_z", &rb_rot_z, 0.0f, 360.0f);

				ImGui::SliderFloat("rb_scale_x", &rb_scale_x, 0.5f, 5.0f);
				ImGui::SliderFloat("rb_scale_y", &rb_scale_y, 0.5f, 5.0f);
				ImGui::SliderFloat("rb_scale_z", &rb_scale_z, 0.5f, 5.0f);

				if (ImGui::Button("ClothJoint"))
				{
					for (int i = 2450; i < 2500; i++) {
						Simulation::m_simul1->getModel()->addRigidBodyParticleBallJoint(0, i);
						Simulation::m_simul2->getModel()->addRigidBodyParticleBallJoint(0, i);
					}
				}

				for (int i = 0; i < 2; i++)
				{
					PBD::Simulation::switchCurrent();
					PBD::SimulationModel* model = PBD::Simulation::getCurrent()->getModel();
					PBD::SimulationModel::RigidBodyVector& rb = model->getRigidBodies();
					if (rb_trans_x != prev_rb_trans_x || rb_trans_y != prev_rb_trans_y || rb_trans_z != prev_rb_trans_z) {

						Vector3r trans_pos = Vector3r(rb_trans_x, rb_trans_y, rb_trans_z);
						rb[n-1]->getLastPosition() = trans_pos;
						rb[n-1]->getOldPosition() = trans_pos;
						rb[n-1]->getPosition() = trans_pos;
						rb[n-1]->getGeometry().updateMeshTransformation(rb[n-1]->getPosition(), rb[n-1]->getRotation().matrix());

						if (i == 1)
						{
							prev_rb_trans_x = rb_trans_x;
							prev_rb_trans_y = rb_trans_y;
							prev_rb_trans_z = rb_trans_z;
						}
					}

					if (rb_rot_x != prev_rb_rot_x)
					{
						rb_rot_w = cos(rb_rot_x / 2);
						rb_rot_x = sin(rb_rot_x / 2);
						rb_rot_y = 0;
						rb_rot_z = 0;

						Quaternionr rot_quaternion = Quaternionr(rb_rot_x, rb_rot_y, rb_rot_z, rb_rot_w);
						rb[n-1]->getLastRotation() = rot_quaternion;
						rb[n-1]->getOldRotation() = rot_quaternion;
						rb[n-1]->getRotation() = rot_quaternion;
						rb[n-1]->getGeometry().updateMeshTransformation(rb[n-1]->getPosition(), rb[n-1]->getRotation().matrix());

						if (i == 1)
						{
							prev_rb_rot_w = rb_rot_w;
							prev_rb_rot_x = rb_rot_x;
							prev_rb_rot_y = rb_rot_y;
							prev_rb_rot_z = rb_rot_z;
						}

					}
					else if (rb_rot_y != prev_rb_rot_y)
					{
						rb_rot_w = cos(rb_rot_y / 2);
						rb_rot_x = 0;
						rb_rot_y = sin(rb_rot_y / 2);
						rb_rot_z = 0;

						Quaternionr rot_quaternion = Quaternionr(rb_rot_x, rb_rot_y, rb_rot_z, rb_rot_w);
						rb[n-1]->getLastRotation() = rot_quaternion;
						rb[n-1]->getOldRotation() = rot_quaternion;
						rb[n-1]->getRotation() = rot_quaternion;
						rb[n-1]->getGeometry().updateMeshTransformation(rb[n-1]->getPosition(), rb[n-1]->getRotation().matrix());

						if (i == 1)
						{
							prev_rb_rot_w = rb_rot_w;
							prev_rb_rot_x = rb_rot_x;
							prev_rb_rot_y = rb_rot_y;
							prev_rb_rot_z = rb_rot_z;
						}
					}
					else if (rb_rot_z != prev_rb_rot_z)
					{
						rb_rot_w = cos(rb_rot_z / 2);
						rb_rot_x = 0;
						rb_rot_y = 0;
						rb_rot_z = sin(rb_rot_z / 2);

						Quaternionr rot_quaternion = Quaternionr(rb_rot_x, rb_rot_y, rb_rot_z, rb_rot_w);
						rb[n-1]->getLastRotation() = rot_quaternion;
						rb[n-1]->getOldRotation() = rot_quaternion;
						rb[n-1]->getRotation() = rot_quaternion;
						rb[n-1]->getGeometry().updateMeshTransformation(rb[n-1]->getPosition(), rb[n-1]->getRotation().matrix());

						if (i == 1)
						{
							prev_rb_rot_w = rb_rot_w;
							prev_rb_rot_x = rb_rot_x;
							prev_rb_rot_y = rb_rot_y;
							prev_rb_rot_z = rb_rot_z;
						}
					}
				}

		


				/*if (rb_scale_x != prev_rb_scale_x || rb_scale_y != prev_rb_scale_y || rb_scale_z != prev_rb_scale_z) {

					rb[n]->getGeometry().updateMeshScale(rb[n]->getPosition(), Vector3r(rb_scale_x / prev_rb_scale_x, rb_scale_y / prev_rb_scale_y, rb_scale_z / prev_rb_scale_z));
					prev_rb_scale_x = rb_scale_x;
					prev_rb_scale_y = rb_scale_y;
					prev_rb_scale_z = rb_scale_z;
				}*/

				//if (ImGui::SliderFloat("x", &m_simulation->m_rigidbodies[n].m_translation(0), -1, 1) ||
				//	ImGui::SliderFloat("y", &m_simulation->m_rigidbodies[n].m_translation(1), -1, 1) ||
				//	ImGui::SliderFloat("z", &m_simulation->m_rigidbodies[n].m_translation(2), -1, 1) ||
				//	ImGui::SliderFloat("rot_angle", &m_simulation->m_rigidbodies[n].m_rotation_angle, -1, 1)
				//	) {
				//	//if (ImGui::IsItemActive())
				//		m_simulation->m_is_simulating = true;
				//}
				/*else
					m_simulation->m_is_simulating = false;*/
				ImGui::EndTabBar();
			}
		}
	}
	ImGui::End();
}

void ImguiManager::createCenterMenu()
{
	ImGui::Begin("center", NULL, window_flags);
	ImGui::SetWindowPos(ImVec2(m_w / 6, 20));
	ImGui::SetWindowSize(ImVec2(m_w / 6 * 3, m_h / 3 * 2), 0);
	if (ImGui::BeginTabBar("center"))
	{
		ImGui::End();
		ImGui::EndTabBar();
	}

}

void ImguiManager::createBottomMenu()
{
	if (ImGui::Begin("Message", NULL, window_flags))
	{
		ImGui::SetWindowPos(ImVec2(180, m_h / 2.5 + 275));
		ImGui::SetWindowSize(ImVec2(2 * m_w / 2.5, 150));
		ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	}
	ImGui::End();

	ImGui::Begin("Files", NULL, window_flags);
	ImGui::SetWindowPos(ImVec2(0, m_h / 2.5 + 425));
	ImGui::SetWindowSize(ImVec2(m_w, 223));
	if (ImGui::BeginTabBar("bottom"))
	{

		if (ImGui::BeginTabItem("Project"))
		{
			createFileDialogBtn();
			if (m_dirPath != "")
			{
				showDirectoryTree(m_dirPath);
			}
		}
	}
	ImGui::End();
}

void ImguiManager::createFileDialogBtn()
{
	if (ImGui::Button("open file dialog"))
		m_fileDialog.Open();


	m_fileDialog.Display();

	if (m_fileDialog.HasSelected())
	{
		std::cout << "Selected filename: " << m_fileDialog.GetSelected().string() << std::endl;
		m_dirPath = m_fileDialog.GetPwd().string();
		std::cout << "Selected filepath: " << m_dirPath << std::endl;

		std::string selected = m_fileDialog.GetSelected().string();
		const size_t last_slash_idx = selected.find_last_of("\\/");
		if (std::string::npos != last_slash_idx)
		{
			selected.erase(0, last_slash_idx + 1);
		}

		if (std::filesystem::path(selected).extension() == ".obj")
		{
			m_filetree_double_clicked_item = selected;

			addRigidbody(m_fileDialog.GetSelected().string().c_str());
			PBD::Simulation::switchCurrent();
			addRigidbody(m_fileDialog.GetSelected().string().c_str());
			PBD::Simulation::switchCurrent();
			//bool res = UtilManager::loadOBJ(m_fileDialog.GetSelected().string().c_str(), vertices, uvs, normals);
		}

		m_fileDialog.ClearSelected();
	}
}

// using file name
void ImguiManager::addRigidbody(std::string fName)
{
	std::string fileName = Utilities::FileSystem::normalizePath(fName);
	Utilities::IndexedFaceMesh mesh;
	PBD::VertexData vd;
	loadObj(fileName, vd, mesh, Vector3r::Ones());

	PBD::SimulationModel* model = PBD::Simulation::getCurrent()->getModel();
	PBD::SimulationModel::RigidBodyVector& rb = model->getRigidBodies();

	PBD::RigidBody* n_rigidbody = new PBD::RigidBody();
	n_rigidbody->initBody(100.0,
		Vector3r(0.0, 5.5, 0.0),
		Quaternionr(1.0, 0.0, 0.0, 0.0),
		vd, mesh,
		Vector3r(5.0, 5.0, 5.0));
	n_rigidbody->setMass(1.0);
	rb.push_back(n_rigidbody);

	const std::vector<Vector3r>* vertices = rb.back()->getGeometry().getVertexDataLocal().getVertices();
	const unsigned int nVert = static_cast<unsigned int>(vertices->size());
	PBD::DistanceFieldCollisionDetection &cd = *(PBD::DistanceFieldCollisionDetection*)PBD::Simulation::getCurrent()->getTimeStep()->getCollisionDetection();
	cd.addCollisionBox(rb.size()-1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector3r(5.0, 5.0, 5.0));

}

// using tree type
void ImguiManager::addRigidbody(treetype tree)
{
	std::string fileName;
	
	switch (tree)
	{
	case cube:
		fileName = Utilities::FileSystem::normalizePath(Utilities::FileSystem::normalizePath(Utilities::FileSystem::getDirectoryPath() + "/" + std::string("data")) + "/models/cube.obj");
		break;
	case cylinder:
		fileName = Utilities::FileSystem::normalizePath(Utilities::FileSystem::normalizePath(Utilities::FileSystem::getDirectoryPath() + "/" + std::string("data")) + "/models/cylinder.obj");
		break;
	case sphere:
		fileName = Utilities::FileSystem::normalizePath(Utilities::FileSystem::normalizePath(Utilities::FileSystem::getDirectoryPath() + "/" + std::string("data")) + "/models/sphere.obj");
		break;
	case cloth:
		break;
	case torus:
		fileName = Utilities::FileSystem::normalizePath(Utilities::FileSystem::normalizePath(Utilities::FileSystem::getDirectoryPath() + "/" + std::string("data")) + "/models/torus.obj");
		break;
	case floor:
		fileName = Utilities::FileSystem::normalizePath(Utilities::FileSystem::normalizePath(Utilities::FileSystem::getDirectoryPath() + "/" + std::string("data")) + "/models/cube.obj");
		break;
	case custom:
		break;
	default:
		break;
	}
	Utilities::IndexedFaceMesh mesh;
	PBD::VertexData vd;
	loadObj(fileName, vd, mesh, Vector3r::Ones());

	PBD::SimulationModel* model = PBD::Simulation::getCurrent()->getModel();
	PBD::SimulationModel::RigidBodyVector& rb = model->getRigidBodies();
	static int s = rb.size();
	PBD::RigidBody* n_rigidbody = new PBD::RigidBody();

	switch (tree)
	{
	case cube:
		n_rigidbody->initBody(1.0,
			Vector3r(5.0, -1.5, 5.0),
			Quaternionr(1.0, 0.0, 0.0, 0.0),
			vd, mesh,
			Vector3r(2.0, 2.0, 2.0));
		n_rigidbody->setMass(0.0);
		n_rigidbody->setFrictionCoeff(static_cast<Real>(0.1));
		break;
	case cylinder:
		n_rigidbody->initBody(1.0,
			Vector3r(5.0, -1.5, 5.0),
			Quaternionr(1.0, 0.0, 0.0, 0.0),
			vd, mesh,
			Vector3r(2.0, 4.0, 2.0));
		n_rigidbody->setMass(1.0);
		n_rigidbody->setFrictionCoeff(static_cast<Real>(0.1));
		break;
	case sphere:
		n_rigidbody->initBody(1.0,
			Vector3r(5.0, -1.5, 5.0),
			Quaternionr(1.0, 0.0, 0.0, 0.0),
			vd, mesh,
			Vector3r(2.0, 2.0, 2.0));
		n_rigidbody->setMass(0.0);
		n_rigidbody->setFrictionCoeff(static_cast<Real>(0.1));
		break;
	case cloth:
		break;
	case torus:
		n_rigidbody->initBody(1.0,
			Vector3r(6.0, -1.5, 5.0),
			Quaternionr(1.0, 0.0, 0.0, 0.0),
			vd, mesh,
			Vector3r(2.0, 2.0, 2.0));
		n_rigidbody->setMass(0.0);
		n_rigidbody->setFrictionCoeff(static_cast<Real>(0.1));
		break;
	case floor:
		n_rigidbody->initBody(1.0,
			Vector3r(0.0, -5.5, 0.0),
			Quaternionr(1.0, 0.0, 0.0, 0.0),
			vd, mesh,
			Vector3r(100.0, 1.0, 100.0));
		n_rigidbody->setMass(0.0);
		break;
	case custom:
		break;
	default:
		break;
	}

	rb.push_back(n_rigidbody);
	//LOG_INFO << "debug" << rb.size() - 1;
	//model->addRigidBodyParticleBallJoint(0, 1);
	const std::vector<Vector3r>* vertices = rb.back()->getGeometry().getVertexDataLocal().getVertices();
	const unsigned int nVert = static_cast<unsigned int>(vertices->size());
	
	PBD::Simulation* current = PBD::Simulation::getCurrent();

	switch (tree)
	{
	case cube:
		if (current == PBD::Simulation::m_simul1) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep()->getCollisionDetection();
			cd.addCollisionBox(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector3r(2.0, 2.0, 2.0));
		}
		else if (current == PBD::Simulation::m_simul2) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep2()->getCollisionDetection();
			cd.addCollisionBox(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector3r(2.0, 2.0, 2.0));
		}
		break;
	case cylinder:
		if (current == PBD::Simulation::m_simul1) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep()->getCollisionDetection();
			cd.addCollisionCylinder(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector2r(2.0, 4.0));
		}
		else if (current == PBD::Simulation::m_simul2) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep2()->getCollisionDetection();
			cd.addCollisionCylinder(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector2r(2.0, 4.0));
		}
		break;
	case sphere:
		if (current == PBD::Simulation::m_simul1) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep()->getCollisionDetection();
			cd.addCollisionSphere(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, 2);
		}
		else if (current == PBD::Simulation::m_simul2) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep2()->getCollisionDetection();
			cd.addCollisionSphere(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, 2);
		}
		break;
	case cloth:
		break;
	case torus:
		if (current == PBD::Simulation::m_simul1) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep()->getCollisionDetection();
			cd.addCollisionTorus(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector2r(2.0, 1.0));
		}
		else if (current == PBD::Simulation::m_simul2) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep2()->getCollisionDetection();
			cd.addCollisionTorus(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector2r(2.0, 1.0));
		}
		break;
	case floor:
		if (current == PBD::Simulation::m_simul1) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep()->getCollisionDetection();
			cd.addCollisionBox(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector3r(100.0, 1.0, 100.0));
		}
		else if (current == PBD::Simulation::m_simul2) {
			PBD::DistanceFieldCollisionDetection& cd = *(PBD::DistanceFieldCollisionDetection*)current->getTimeStep2()->getCollisionDetection();
			cd.addCollisionBox(rb.size() - 1, PBD::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType, &(*vertices)[0], nVert, Vector3r(100.0, 1.0, 100.0));
		}
		break;
	case custom:
		break;
	default:
		break;
	}
}

void ImguiManager::ItemRowsBackground(float lineHeight, const ImColor& color)
{
	auto* drawList = ImGui::GetWindowDrawList();
	const auto& style = ImGui::GetStyle();

	if (lineHeight < 0)
	{
		lineHeight = ImGui::GetTextLineHeight();
	}
	lineHeight += style.ItemSpacing.y;

	float scrollOffsetH = ImGui::GetScrollX();
	float scrollOffsetV = ImGui::GetScrollY();
	float scrolledOutLines = floorf(scrollOffsetV / lineHeight);
	scrollOffsetV -= lineHeight * scrolledOutLines;

	ImVec2 clipRectMin(ImGui::GetWindowPos().x, ImGui::GetWindowPos().y);
	ImVec2 clipRectMax(clipRectMin.x + ImGui::GetWindowWidth(), clipRectMin.y + ImGui::GetWindowHeight());

	if (ImGui::GetScrollMaxX() > 0)
	{
		clipRectMax.y -= style.ScrollbarSize;
	}

	drawList->PushClipRect(clipRectMin, clipRectMax);

	bool isOdd = (static_cast<int>(scrolledOutLines) % 2) == 0;

	float yMin = clipRectMin.y - scrollOffsetV + ImGui::GetCursorPosY();
	float yMax = clipRectMax.y - scrollOffsetV + lineHeight;
	float xMin = clipRectMin.x + scrollOffsetH + ImGui::GetWindowContentRegionMin().x;
	float xMax = clipRectMin.x + scrollOffsetH + ImGui::GetWindowContentRegionMax().x;

	for (float y = yMin; y < yMax; y += lineHeight, isOdd = !isOdd)
	{
		if (isOdd)
		{
			drawList->AddRectFilled({ xMin, y - style.ItemSpacing.y }, { xMax, y + lineHeight }, color);
		}
	}

	drawList->PopClipRect();
}

std::pair<bool, uint32_t> ImguiManager::DirectoryTreeViewRecursive(const std::filesystem::path& path, uint32_t* count, int* selection_mask)
{
	ImGuiTreeNodeFlags base_flags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick | ImGuiTreeNodeFlags_SpanAvailWidth | ImGuiTreeNodeFlags_SpanFullWidth;

	bool any_node_clicked = false;
	uint32_t node_clicked = 0;

	for (const auto& entry : std::filesystem::directory_iterator(path))
	{
		ImGuiTreeNodeFlags node_flags = base_flags;
		const bool is_selected = (*selection_mask & BIT(*count)) != 0;
		if (is_selected)
			node_flags |= ImGuiTreeNodeFlags_Selected;

		std::string name = entry.path().string();

		auto lastSlash = name.find_last_of("/\\");
		lastSlash = lastSlash == std::string::npos ? 0 : lastSlash + 1;
		name = name.substr(lastSlash, name.size() - lastSlash);

		bool entryIsFile = !std::filesystem::is_directory(entry.path());
		if (entryIsFile)
			node_flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;

		bool node_open = ImGui::TreeNodeEx((void*)(intptr_t)(*count), node_flags, name.c_str());

		if (ImGui::IsItemClicked())
		{
			node_clicked = *count;
			any_node_clicked = true;
			if (ImGui::IsMouseDoubleClicked(0) && ImGui::IsItemHovered() && std::filesystem::path(name).extension() == ".obj")
			{
				std::cout << name << std::endl;
				m_filetree_double_clicked_item = name.c_str();

				addRigidbody(entry.path().string());
				PBD::Simulation::switchCurrent();
				addRigidbody(entry.path().string());
				PBD::Simulation::switchCurrent();
				/*std::string fileName = Utilities::FileSystem::normalizePath(entry.path().string().c_str());
				Utilities::IndexedFaceMesh mesh;
				PBD::VertexData vd;
				loadObj(fileName, vd, mesh, Vector3r::Ones());*/

				//bool res = UtilManager::loadOBJ(entry.path().string().c_str(), vertices, uvs, normals);
			}
		}

		(*count)--;

		if (!entryIsFile)
		{
			if (node_open)
			{

				auto clickState = DirectoryTreeViewRecursive(entry.path(), count, selection_mask);

				if (!any_node_clicked)
				{
					any_node_clicked = clickState.first;
					node_clicked = clickState.second;
				}

				ImGui::TreePop();
			}
			else
			{
				if (ImGui::IsMouseDragging && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem))
				{
					//m_MovePath = dirInfo.absolutePath.c_str();
				}
				for (const auto& e : std::filesystem::recursive_directory_iterator(entry.path()))
					(*count)--;
			}
		}
	}

	return { any_node_clicked, node_clicked };
}

void ImguiManager::showDirectoryTree(std::string directoryPath)
{
	ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2{ 0.0f, 0.0f });

	{
		uint32_t count = 0;


		for (const auto& entry : std::filesystem::recursive_directory_iterator(directoryPath)) {
			count++;
		}

		static int selection_mask = 0;

		auto clickState = DirectoryTreeViewRecursive(directoryPath, &count, &selection_mask);

		if (clickState.first)
		{
			// Update selection state
			if (ImGui::GetIO().KeyCtrl)
				selection_mask ^= BIT(clickState.second);
			else
				selection_mask = BIT(clickState.second);
		}
	}

	ImGui::PopStyleVar();
}

void ImguiManager::fbo_init()
{
	int window_width = m_w;
	int window_height = m_h;

	float rotation_degree = 0.0f;

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);

	glGenFramebuffersEXT(1, &m_fbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fbo); // Bind our frame buffer  

	glGenRenderbuffersEXT(1, &m_fbo_depth);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, m_fbo_depth);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, window_width, window_height);
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_fbo_depth);

	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);

	glGenTextures(1, &m_fbo_texture);
	glBindTexture(GL_TEXTURE_2D, m_fbo_texture);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, window_width, window_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Unbind the texture  
	glBindTexture(GL_TEXTURE_2D, 0);

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_fbo_texture, 0); // Attach the texture fbo_texture to the color buffer in our frame buffer  

	//glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_fbo_depth); // Attach the depth buffer fbo_depth to our frame buffer  

	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);

	if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
	{
		std::cout << "Couldn't create frame buffer" << std::endl;
		exit(0);
	}
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0); // Unbind our frame buffer 
}

void ImguiManager::fbo2_init()
{
	int window_width = m_w;
	int window_height = m_h;

	float rotation_degree = 0.0f;

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_DEPTH_TEST);

	glGenFramebuffersEXT(1, &m_fbo2);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fbo2); // Bind our frame buffer  

	glGenRenderbuffersEXT(1, &m_fbo_depth2);
	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, m_fbo_depth2);
	glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, window_width, window_height);
	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_fbo_depth2);

	glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);



	glGenTextures(1, &m_fbo_texture2);
	glBindTexture(GL_TEXTURE_2D, m_fbo_texture2);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, window_width, window_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Unbind the texture  
	glBindTexture(GL_TEXTURE_2D, 0);


	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_fbo_texture2, 0); // Attach the texture fbo_texture to the color buffer in our frame buffer  


	//glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_fbo_depth); // Attach the depth buffer fbo_depth to our frame buffer  

	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);

	if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
	{
		std::cout << "Couldn't create frame buffer" << std::endl;
		exit(0);
	}
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0); // Unbind our frame buffer 
}

void ImguiManager::fbo_bind(unsigned int n)
{
	switch (n)
	{
	case 0:
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fbo); // Bind our frame buffer for rendering
		break;
	case 1:
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fbo2); // Bind our frame buffer for rendering 
		break;
	default:
		break;
	}
	

	glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
	glViewport(0, 0, m_w, m_h);


	glClearColor(0.53f, 0.51f, 0.51f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();


}
void ImguiManager::fbo_unbind(unsigned int n)
{
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0); // Unbind our texture 

	glClearColor(1.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	switch (n)
	{
	case 0:
		glBindTexture(GL_TEXTURE_2D, m_fbo_texture); // Bind our frame buffer texture 
		break;
	case 1:
		glBindTexture(GL_TEXTURE_2D, m_fbo_texture2); // Bind our frame buffer texture 
		break;
	default:
		break;
	}
}

void ImguiManager::fbo_cleanup()
{
	glDeleteTextures(1, &m_fbo_texture);
	glDeleteRenderbuffersEXT(1, &m_fbo_depth);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	glDeleteFramebuffersEXT(1, &m_fbo);

	glDeleteTextures(1, &m_fbo_texture2);
	glDeleteRenderbuffersEXT(1, &m_fbo_depth2);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	glDeleteFramebuffersEXT(1, &m_fbo2);
}


void ImguiManager::loadObj(const std::string& filename, PBD::VertexData& vd, Utilities::IndexedFaceMesh& mesh, const Vector3r& scale)
{
	std::vector<Utilities::OBJLoader::Vec3f> x;
	std::vector<Utilities::OBJLoader::Vec3f> normals;
	std::vector<Utilities::OBJLoader::Vec2f> texCoords;
	std::vector<Utilities::MeshFaceIndices> faces;
	Utilities::OBJLoader::Vec3f s = { (float)scale[0], (float)scale[1], (float)scale[2] };
	Utilities::OBJLoader::loadObj(filename, &x, &faces, &normals, &texCoords, s);

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
