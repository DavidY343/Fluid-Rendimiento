//
// Created by david on 11/21/23.
//
#include "gtest/gtest.h"
#include <cstdlib>  // Para la función std::system
#include <filesystem>
class TestFuncionales : public ::testing::Test{};

//comprueba que la trza small-1.fld es la misma que sale cuando le paso small.fld y 1 iteracion
TEST_F(TestFuncionales, testSmall1) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/small.fld";
  std::string const ruta_salida = "../../out/small-1.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/small-1.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(1) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza small-2.fld es la misma que sale cuando le paso small.fld y 2 iteracion
TEST_F(TestFuncionales, testSmall2) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/small.fld";
  std::string const ruta_salida = "../../out/small-2.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/small-2.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(2) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza small-3.fld es la misma que sale cuando le paso small.fld y 3 iteraciones
TEST_F(TestFuncionales, testSmall3) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/small.fld";
  std::string const ruta_salida = "../../out/small-3.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/small-3.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(3) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza small-4.fld es la misma que sale cuando le paso small.fld y 4 iteraciones
TEST_F(TestFuncionales, testSmall4) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/small.fld";
  std::string const ruta_salida = "../../out/small-4.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/small-4.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(4) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza small-5.fld es la misma que sale cuando le paso small.fld y 5 iteraciones
TEST_F(TestFuncionales, testSmall5) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/small.fld";
  std::string const ruta_salida = "../../out/small-5.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/small-5.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(5) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}
//comprueba que la trza large-1.fld es la misma que sale cuando le paso large.fld y 1 iteraciones
TEST_F(TestFuncionales, testlarge1) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/large.fld";
  std::string const ruta_salida = "../../out/large-1.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/large-1.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(1) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza large-2.fld es la misma que sale cuando le paso large.fld y 2 iteraciones
TEST_F(TestFuncionales, testlarge2) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/large.fld";
  std::string const ruta_salida = "../../out/large-2.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/large-2.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(2) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza large-3.fld es la misma que sale cuando le paso large.fld y 3 iteraciones
TEST_F(TestFuncionales, testlarge3) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/large.fld";
  std::string const ruta_salida = "../../out/large-3.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/large-3.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(3) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza large-4.fld es la misma que sale cuando le paso large.fld y 4 iteraciones
TEST_F(TestFuncionales, testlarge4) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/large.fld";
  std::string const ruta_salida = "../../out/large-4.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/large-4.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(4) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}

//comprueba que la trza large-5.fld es la misma que sale cuando le paso large.fld y 5 iteraciones
TEST_F(TestFuncionales, testlarge5) {
  // Definir las rutas de los archivos de prueba y salida
  std::string const ruta_entrada = "../../in/large.fld";
  std::string const ruta_salida = "../../out/large-5.fld";
  std::filesystem::path ruta_ejecutable = std::filesystem::current_path();  // Ruta del directorio actual
  ruta_ejecutable /= "../fluid/fluid";  // Agregar "fluid" al final de la ruta
  std::string const ruta_traza = "../../trz/large-5.fld";
  // Construir el comando para ejecutar el programa principal con rutas absolutas
  std::string const comando = ruta_ejecutable.string() + " " + std::to_string(5) + " " + ruta_entrada + " " + ruta_salida;

  // Ejecutar el comando
  int const resultado = std::system(comando.c_str());

  // Comprobar que la ejecución fue exitosa (código de salida 0)
  EXPECT_EQ(resultado, 0);
  std::string const diff = std::string("diff") + " " + ruta_salida + " " + ruta_traza;
  int const final = std::system(diff.c_str());
  EXPECT_EQ(final, 0);
}