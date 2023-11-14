//
// Created by david on 11/7/23.
//

#include "progargs.hpp"

#include "grid.hpp"
//
// Created by david on 10/3/23.
//

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Funcion que comprueba que los archivos se abrieron correctamente
int comprobar_archivos(std::vector<std::string> const & argumentos, std::ifstream const & inputFile,
                       std::ofstream const & outputFile) {
  if (!inputFile.is_open()) {  // Si nos devuelve False, error
    std::cerr << "Error: Cannot open " << argumentos[1] << " for reading" << '\n';
    return -3;
  }
  if (!outputFile.is_open()) {  // Si nos devuelve False, error
    std::cerr << "Error: Cannot open " << argumentos[2] << " for writing" << '\n';
    return -4;
  }
  return 0;
}

// Si iteraciones es un numero
int comprobar_params(std::vector<std::string> const & argumentos,
                     std::ifstream const & inputFile, std::ofstream const & outputFile) {
  int nSteps = 0;
  try {
    nSteps = std::stoi(argumentos[0]);
  } catch (
      std::invalid_argument const & excepcion) {  // Solo ocurrira en caso de que argv[1] no sea int
    std::cerr << "Error: timesteps must be numeric."
              << "\n";
    return -1;
  }
  if (nSteps < 0) {
    std::cerr << "Error: Invalid number of time steps."
              << "\n";
    return -2;
  }
  return comprobar_archivos(argumentos, inputFile, outputFile);
};

int n_params(int argc){
  if (argc != 4) {
    std::cerr << "Error: Invalid number of arguments:" << argc << "\n";
    return -1;
  }
  return 0;
}

// Funciones para leer del archivo
template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char * as_writable_buffer(T & value) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  return reinterpret_cast<char *>(&value);
}

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char const * as_buffer(T const & value) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  return reinterpret_cast<char const *>(&value);
}

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
T read_binary_value(std::istream & is) {
  T value{};
  is.read(as_writable_buffer(value), sizeof(value));
  return value;
}

void comprobar_fallos_cabecera(std::vector<particula> const & particulas, int n_particulas_int) {
  int const longitud = static_cast<int>(particulas.size());
  int const error    = -5;
  if (n_particulas_int <= 0) {
    std::cerr << "Invalid number of particles: " << n_particulas_int<< ".\n";
    exit(error);
  }
  if (longitud != n_particulas_int) {
    std::cerr << "Number of particles mismatch. Header: " << n_particulas_int
              << ", Found: " << longitud << ".\n";
    exit(error);
  }
}

// Con esta funcion creamos el array de partículas, y generamos los datos del archivo de entrada
std::vector<particula> crear_particulas(std::ifstream const & inputFile) {
  // Obtenemos el número de partículas de la cabecera
  auto n_particulas_int = read_binary_value<int>((std::istream &) inputFile);
  std::vector<particula> particulas;  // Creamos el vector con las partículas
  int ident = 0;
  while (!inputFile.eof()) {
    // Leer los datos de una partícula
    int const id_dat = ident;
    auto px_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto py_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto pz_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto hvx         = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto hvy         = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto hvz         = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto vx_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto vy_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto vz_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    if (inputFile.eof()) { break; }
    particula const n_particula(id_dat, px_dat, py_dat, pz_dat, hvx, hvy, hvz, vy_dat, vx_dat,
                                vz_dat);
    particulas.push_back(n_particula);
    ident++;
  }
  comprobar_fallos_cabecera(particulas, n_particulas_int);
  return particulas;
}

void particula::colisionLimiteEjeX(bool lim_inf) {
  double const min_value = 0.0000000001;
  if (lim_inf) {
    double const difLimX = constantes::dp_const - (getpx() - constantes::bmin_const[0]);
    if (difLimX > min_value) {
      setax(getax() +
            (constantes::sc_const * constantes::t_const - constantes::dv_const * getvx()));
    }
  } else {
    double const difLimX = constantes::dp_const - (constantes::bmax_const[0] - getpx());
    if (difLimX > min_value) {
      setax(getax() - (constantes::sc_const * difLimX + constantes::dv_const * getvx()));
    }
  }
}

void particula::colisionLimiteEjeY(bool lim_inf) {
  double const min_value = 0.0000000001;
  if (lim_inf) {
    double const difLimY = constantes::dp_const - (getpy() - constantes::bmin_const[1]);
    if (difLimY > min_value) {
      setay(getay() +
            (constantes::sc_const * constantes::t_const - constantes::dv_const * getvy()));
    }
  } else {
    double const difLimY = constantes::dp_const - (constantes::bmax_const[1] - getpy());
    if (difLimY > min_value) {
      setay(getay() - (constantes::sc_const * difLimY + constantes::dv_const * getvy()));
    }
  }
}

void particula::colisionLimiteEjeZ(bool lim_inf) {
  double const min_value = 0.0000000001;
  if (lim_inf) {
    double const difLimZ = constantes::dp_const - (getpz() - constantes::bmin_const[2]);
    if (difLimZ > min_value) {
      setaz(getaz() +
            (constantes::sc_const * constantes::t_const - constantes::dv_const * getvz()));
    }
  } else {
    double const difLimZ = constantes::dp_const - (constantes::bmax_const[2] - getpz());
    if (difLimZ > min_value) {
      setaz(getaz() - (constantes::sc_const * difLimZ + constantes::dv_const * getvz()));
    }
  }
}

void particula::limiteRecintox(bool lim_inf) {
  double const dx_dat =
      lim_inf ? (getpx() - constantes::bmin_const[0]) : (constantes::bmax_const[0] - getpx());
  if (dx_dat < 0) {
    px  = lim_inf ? (constantes::bmin_const[0] - dx_dat) : (constantes::bmax_const[0] + dx_dat);
    vx  = -vx;
    hvx = -hvx;
  }
}

void particula::limiteRecintoy(bool lim_inf) {
  double const dy_dat =
      lim_inf ? (getpy() - constantes::bmin_const[1]) : (constantes::bmax_const[1] - getpy());
  if (dy_dat < 0) {
    py  = lim_inf ? (constantes::bmin_const[1] - dy_dat) : (constantes::bmax_const[1] + dy_dat);
    vy  = -vy;
    hvy = -hvy;
  }
}

void particula::limiteRecintoz(bool lim_inf) {
  double const dz_dat =
      lim_inf ? (getpz() - constantes::bmin_const[2]) : (constantes::bmax_const[2] - getpz());
  if (dz_dat < 0) {
    pz  = lim_inf ? (constantes::bmin_const[2] - dz_dat) : (constantes::bmax_const[2] + dz_dat);
    vz  = -vz;
    hvz = -hvz;
  }
}

void particula::actualizarMovimiento() {
  setpx(getpx() + gethvx() * constantes::t_const +
        getax() * (constantes::t_const * constantes::t_const));
  setvx(gethvx() + (getax() * constantes::t_const) / 2);
  sethvx(gethvx() + getax() * constantes::t_const);

  setpy(getpy() + gethvy() * constantes::t_const +
        getay() * (constantes::t_const * constantes::t_const));
  setvy(gethvy() + (getay() * constantes::t_const) / 2);
  sethvy(gethvy() + getay() * constantes::t_const);

  setpz(getpz() + gethvz() * constantes::t_const +
        getaz() * (constantes::t_const * constantes::t_const));
  setvz(gethvz() + (getaz() * constantes::t_const) / 2);
  sethvz(gethvz() + getaz() * constantes::t_const);
}