# demo 目录: 以 ctest 验收用例为主 (每个子目录是一个独立测试程序)。
# apps/ 是附加区 —— 现象演示, 不注册 ctest (见 apps/CMakeLists.txt)。
add_subdirectory(core)
add_subdirectory(autodiff)
add_subdirectory(linalg)
add_subdirectory(apps)
