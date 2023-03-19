from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout


class helloRecipe(ConanFile):
    name = "mct-dynamics"
    version = "0.0"

    # Optional metadata
    author = "Jonas Nußdorfer jonas.nussdorfer@gmail.com"
    url = "https://github.com/nussjo/mct-dynamics>"
    description = "Implementation of solving the MCT dynamics equations"
    topics = ("MCT")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "src/*", "include/*"

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def requirements(self):
        self.requires('etl/20.35.11@')

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["hello"]