from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.scm import Git


class MctDynamics(ConanFile):
    name = "mct-dynamics"
    version = "0.0"

    # Optional metadata
    author = "Jonas Nußdorfer jonas.nussdorfer@gmail.com"
    url = "https://github.com/nussjo/mct-dynamics"
    description = "Implementation of solving the MCT dynamics equations"
    license = "GPL-3.0"
    topics = ("MCT")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def source(self):
        git = Git(self)
        git.clone(url=self.url, target=".")
        # git.checkout(self.version)

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
        pass