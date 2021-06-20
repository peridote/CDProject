#include "Utils/OBJLoader.h"
#include "Simulation/SimulationModel.h"
#include "Utils/FileSystem.h"
#include "Simulation/Simulation.h"
#include "math.h"


#include "ImguiManager.h"


ImguiManager::ImguiManager()
{
	m_window = NULL;
	m_w = m_h = 0;
	m_clear_color = ImVec4(0.0f, 0.0f, 0.0f, 0.0f);
	m_filetreeitem_current_idx = 0;
	m_filetree_double_clicked_item = "";
	m_filetreeitems = {"(1)torus.obj"};
	m_filetree_num = 2; // cloth collision demo 에 이미 2개 생성 되있음

	translation = Eigen::Vector3f(0, 0, 0);
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
	m_clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

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

	fbo_init();
}

void ImguiManager::StartFrame()
{
	glfwPollEvents();

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
	m_filetreeitems.clear();
	m_filetreeitems.push_back("(1)torus.obj"); // 현재 collision demo 기본 obj가 tours라서 넣어줌...
	m_filetree_num = 2;
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

void ImguiManager::createLeftSideMenu()
{
	
	if (ImGui::Begin("DockSpace"))
	{
		//ImGui::SetWindowSize(ImVec2(m_w / 6, m_h / 3 * 2));
		ImGui::BeginDockspace();

		// left Menu
		if (ImGui::BeginDock("Project Files"))
		{
			//ImGui::SetWindowPos(ImVec2(0, 20));
			//ImGui::SetWindowSize(ImVec2(m_w / 6, m_h / 3 * 2));

			if (m_filetree_double_clicked_item != "")
			{
				m_filetree_double_clicked_item = "(" + std::to_string(m_filetree_num) + ")" + m_filetree_double_clicked_item;
				m_map_filetree_rigidbody.insert(std::make_pair(m_filetree_double_clicked_item, m_filetree_num++));
				m_filetreeitems.push_back(m_filetree_double_clicked_item);
				std::cout << m_filetree_double_clicked_item << std::endl;
				m_filetree_double_clicked_item = "";
			}

			for (int n = 0; n < m_filetreeitems.size(); n++)
			{
				ImGui::PushID(n);
				const bool is_selected = (m_filetreeitem_current_idx == n);
				if (ImGui::Selectable(m_filetreeitems[n].c_str(), is_selected)) {
					if (n != m_filetreeitem_current_idx)
					{
						m_filetree_current_item = m_filetreeitems[n];
						m_filetreeitem_current_idx = n;
					}
				}

				if (is_selected)
					ImGui::SetItemDefaultFocus();
				ImGui::PopID();
			}
		}
		ImGui::EndDock();

		// rightMenu
		if (ImGui::BeginDock("Properties"))
		{
			//ImGui::SetWindowPos(ImVec2(m_w / 6 * 4, 20));
			//ImGui::SetWindowSize(ImVec2(m_w / 6 * 2, m_h / 3 * 2), 0);
			
			if (ImGui::BeginTabBar("right"))
			{
				
				if (strcmp(m_filetree_current_item.c_str(), "(1)torus.obj") == 0)
				{
					static bool test1 = true;
					ImGui::Checkbox("Collision", &test1);
					static bool test2 = true;
					ImGui::Checkbox("Visible", &test2);
					static bool test3 = true;
					ImGui::Checkbox("Static", &test3);
					ImGui::EndTabBar();
				}

				else if (std::filesystem::path(m_filetree_current_item).extension() == ".obj")
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

					static float rb_trans_x = rb[n]->getPosition()(0);
					static float rb_trans_y = rb[n]->getPosition()(1);
					static float rb_trans_z = rb[n]->getPosition()(2);
					static float prev_rb_trans_x = 0.0f;
					static float prev_rb_trans_y = 0.0f;
					static float prev_rb_trans_z = 0.0f;

					static float rb_rot_w = rb[n]->getRotation().w();
					static float rb_rot_x = rb[n]->getRotation().x();
					static float rb_rot_y = rb[n]->getRotation().y();
					static float rb_rot_z = rb[n]->getRotation().z();
					static float prev_rb_rot_w = 1.0f;
					static float prev_rb_rot_x = 0.0f;
					static float prev_rb_rot_y = 0.0f;
					static float prev_rb_rot_z = 0.0f;


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

					if (rb_trans_x != prev_rb_trans_x || rb_trans_y != prev_rb_trans_y || rb_trans_z != prev_rb_trans_z) {
						
						Vector3r trans_pos = Vector3r(rb_trans_x, rb_trans_y, rb_trans_z);
						rb[n]->getLastPosition() = trans_pos;
						rb[n]->getOldPosition() = trans_pos;
						rb[n]->getPosition() = trans_pos;
						rb[n]->getGeometry().updateMeshTransformation(rb[n]->getPosition(), rb[n]->getRotation().matrix());

						prev_rb_trans_x = rb_trans_x;
						prev_rb_trans_y = rb_trans_y;
						prev_rb_trans_z = rb_trans_z;

						
					}

					if (rb_rot_x != prev_rb_rot_x)
					{
						rb_rot_w = cos(rb_rot_x / 2);
						rb_rot_x = sin(rb_rot_x/2);
						rb_rot_y = 0;
						rb_rot_z = 0;

						Quaternionr rot_quaternion = Quaternionr(rb_rot_x, rb_rot_y, rb_rot_z, rb_rot_w);
						rb[n]->getLastRotation() = rot_quaternion;
						rb[n]->getOldRotation() = rot_quaternion;
						rb[n]->getRotation() = rot_quaternion;
						rb[n]->getGeometry().updateMeshTransformation(rb[n]->getPosition(), rb[n]->getRotation().matrix());

						prev_rb_rot_w = rb_rot_w;
						prev_rb_rot_x = rb_rot_x;
						prev_rb_rot_y = rb_rot_y;
						prev_rb_rot_z = rb_rot_z;
					}
					else if (rb_rot_y != prev_rb_rot_y)
					{ 
						rb_rot_w = cos(rb_rot_y / 2);
						rb_rot_x = 0;
						rb_rot_y = sin(rb_rot_y / 2);
						rb_rot_z = 0;

						Quaternionr rot_quaternion = Quaternionr(rb_rot_x, rb_rot_y, rb_rot_z, rb_rot_w);
						rb[n]->getLastRotation() = rot_quaternion;
						rb[n]->getOldRotation() = rot_quaternion;
						rb[n]->getRotation() = rot_quaternion;
						rb[n]->getGeometry().updateMeshTransformation(rb[n]->getPosition(), rb[n]->getRotation().matrix());
						prev_rb_rot_w = rb_rot_w;
						prev_rb_rot_x = rb_rot_x;
						prev_rb_rot_y = rb_rot_y;
						prev_rb_rot_z = rb_rot_z;
					}
					else if (rb_rot_z != prev_rb_rot_z)
					{ 
						rb_rot_w = cos(rb_rot_z / 2);
						rb_rot_x = 0;
						rb_rot_y = 0;
						rb_rot_z = sin(rb_rot_z / 2);

						Quaternionr rot_quaternion = Quaternionr(rb_rot_x, rb_rot_y, rb_rot_z, rb_rot_w);
						rb[n]->getLastRotation() = rot_quaternion;
						rb[n]->getOldRotation() = rot_quaternion;
						rb[n]->getRotation() = rot_quaternion;
						rb[n]->getGeometry().updateMeshTransformation(rb[n]->getPosition(), rb[n]->getRotation().matrix());

						prev_rb_rot_w = rb_rot_w;
						prev_rb_rot_x = rb_rot_x;
						prev_rb_rot_y = rb_rot_y;
						prev_rb_rot_z = rb_rot_z;
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
		ImGui::EndDock();
		
		ImGui::EndDockspace();
	}
	ImGui::End();
}

void ImguiManager::createRightSideMenu()
{
	if (ImGui::Begin("right"))
	{
		ImGui::BeginDockspace();
		if (ImGui::BeginDock("right"))
		{
			//ImGui::SetWindowPos(ImVec2(m_w / 6 * 4, 20));
			//ImGui::SetWindowSize(ImVec2(m_w / 6 * 2, m_h / 3 * 2), 0);
			if (ImGui::BeginTabBar("right"))
			{
				if (strcmp(m_filetree_current_item.c_str(), "a") == 0)
				{
					static bool test1 = true;
					ImGui::Checkbox("Test", &test1);
					static bool test2 = true;
					ImGui::Checkbox("Test2", &test2);
					ImGui::EndTabBar();
				}

				else if (std::filesystem::path(m_filetree_current_item).extension() == ".obj")
				{
					int n = m_map_filetree_rigidbody.find(m_filetree_current_item)->second;
					static bool test1 = true;
					ImGui::Checkbox("Test", &test1);
					static bool test2 = true;
					ImGui::Checkbox("Test2", &test2);
					static bool test3 = true;
					ImGui::Checkbox("Test3", &test3);
					static bool test4 = true;
					ImGui::Checkbox("Test4", &test4);

					PBD::SimulationModel* model = PBD::Simulation::getCurrent()->getModel();
					PBD::SimulationModel::RigidBodyVector& rb = model->getRigidBodies();

					static float rb_scale_x = 1.0f;
					static float rb_scale_y = 1.0f;
					static float rb_scale_z = 1.0f;
					static float prev_rb_scale_x = 1.0f;
					static float prev_rb_scale_y = 1.0f;
					static float prev_rb_scale_z = 1.0f;

					static float rb_trans_x = rb[n]->getPosition()(0);
					static float rb_trans_y = rb[n]->getPosition()(1);
					static float rb_trans_z = rb[n]->getPosition()(2);
					static float prev_rb_trans_x = 0.0f;
					static float prev_rb_trans_y = 0.0f;
					static float prev_rb_trans_z = 0.0f;

					static float rb_rot_w = rb[n]->getRotation().w();
					static float rb_rot_x = rb[n]->getRotation().x();
					static float rb_rot_y = rb[n]->getRotation().y();
					static float rb_rot_z = rb[n]->getRotation().z();
					static float prev_rb_rot_w = 0.0f;
					static float prev_rb_rot_x = 0.0f;
					static float prev_rb_rot_y = 0.0f;
					static float prev_rb_rot_z = 0.0f;

					ImGui::SliderFloat("rb_trans_x", &rb_trans_x, -5.0f, 5.0f);
					ImGui::SliderFloat("rb_trans_y", &rb_trans_y, -5.0f, 5.0f);
					ImGui::SliderFloat("rb_trans_z", &rb_trans_z, -5.0f, 5.0f);

					ImGui::SliderFloat("rb_rot_x", &rb_rot_w, -1.0f, 1.0f);
					ImGui::SliderFloat("rb_rot_x", &rb_rot_x, -1.0f, 1.0f);
					ImGui::SliderFloat("rb_rot_y", &rb_rot_y, -1.0f, 1.0f);
					ImGui::SliderFloat("rb_rot_z", &rb_rot_z, -1.0f, 1.0f);

					ImGui::SliderFloat("rb_scale_x", &rb_scale_x, 0.5f, 5.0f);
					ImGui::SliderFloat("rb_scale_y", &rb_scale_y, 0.5f, 5.0f);
					ImGui::SliderFloat("rb_scale_z", &rb_scale_z, 0.5f, 5.0f);

					//if (rb_trans_x != prev_rb_trans_x || rb_trans_y != prev_rb_trans_y || rb_trans_z != prev_rb_trans_z
					//	|| rb_rot_x != prev_rb_rot_x || rb_rot_y != prev_rb_rot_y || rb_rot_z != prev_rb_rot_z) {

					//	rb[n]->getGeometry().updateMeshTransformation(rb[n]->getPosition() + Vector3r(rb_trans_x, rb_trans_y, rb_trans_z), Quaternionr(rb_rot_x, rb_rot_y, rb_rot_z, rb_rot_w).matrix());

					//	prev_rb_trans_x = rb_trans_x;
					//	prev_rb_trans_y = rb_trans_y;
					//	prev_rb_trans_z = rb_trans_z;

					//	//prev_rb_rot_w = rb_rot_w;
					//	prev_rb_rot_x = rb_rot_x;
					//	prev_rb_rot_y = rb_rot_y;
					//	prev_rb_rot_z = rb_rot_z;
					//}

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
		ImGui::EndDock();
		ImGui::EndDockspace();
	}
	ImGui::End();

}

void ImguiManager::createCenterMenu()
{
	ImGui::Begin("center");
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
	ImGui::Begin("Files");
	ImGui::SetWindowPos(ImVec2(0, m_h / 3 * 2 + 140));
	ImGui::SetWindowSize(ImVec2(m_w, m_h / 3 - 180));
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
			//bool res = UtilManager::loadOBJ(m_fileDialog.GetSelected().string().c_str(), vertices, uvs, normals);
		}

		m_fileDialog.ClearSelected();
	}
}

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
	//glEnable(GL_DEPTH_TEST);
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

	glGenFramebuffersEXT(1, &m_fbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fbo); // Bind our frame buffer  

	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, m_fbo_texture, 0); // Attach the texture fbo_texture to the color buffer in our frame buffer  

	glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, m_fbo_depth); // Attach the depth buffer fbo_depth to our frame buffer  

	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);

	if (status != GL_FRAMEBUFFER_COMPLETE_EXT)
	{
		std::cout << "Couldn't create frame buffer" << std::endl;
		exit(0);
	}
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0); // Unbind our frame buffer 
}
void ImguiManager::fbo_bind()
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_fbo); // Bind our frame buffer for rendering 

	glPushAttrib(GL_VIEWPORT_BIT | GL_ENABLE_BIT);
	glViewport(0, 0, m_w, m_h);


	glClearColor(0.53f, 0.51f, 0.51f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();


}
void ImguiManager::fbo_unbind()
{
	glPopAttrib();
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0); // Unbind our texture 

	glClearColor(1.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();


	glBindTexture(GL_TEXTURE_2D, m_fbo_texture); // Bind our frame buffer texture 
}

void ImguiManager::fbo_cleanup()
{
	glDeleteTextures(1, &m_fbo_texture);
	glDeleteRenderbuffersEXT(1, &m_fbo_depth);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	glDeleteFramebuffersEXT(1, &m_fbo);
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
