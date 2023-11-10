//
// Created by david on 11/7/23.
//

#include "progargs.hpp"
#include "grid.hpp"
//
// Created by david on 10/3/23.
//

#include <cmath>
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

// Funcion que comprueba el número de parámetros pasados y si iteraciones es un numero
int comprobar_params(int argc, std::vector<std::string> const & argumentos,
                     std::ifstream const & inputFile, std::ofstream const & outputFile) {
  if (argc != 4) {
    std::cerr << "Error: Invalid number of arguments:" << argc << "\n";
    return -1;
  }
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
    std::cerr << "Error: Invalid number of timesteps."
              << "\n";
    return -2;
  }
  return comprobar_archivos(argumentos, inputFile, outputFile);
};

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
/*
std::vector<double> calculo_tmñ_malla(std::vector<double> const & max_min, std::vector<double> const & max, double h) {
  double const nx_vect = (max_min[0] - max_min[1]) / h;
  double const ny_vect = (max_min[2] - max_min[3]) / h;
  double const nz_vect = (max_min[4] - max_min[5]) / h;

  return {nx_vect, ny_vect, nz_vect};
}
*/
/*
std::vector<double> calcular_max_min(std::vector<particula> const & particulas) {
  double max_px = -std::numeric_limits<double>::max();
  double min_px = std::numeric_limits<double>::max();
  double max_py = -std::numeric_limits<double>::max();
  double min_py = std::numeric_limits<double>::max();
  double max_pz = -std::numeric_limits<double>::max();
  double min_pz = std::numeric_limits<double>::max();

  for (particula const & n_particula : particulas) {
    double const px_dat = n_particula.getpx();
    double const py_dat = n_particula.getpy();
    double const pz_dat = n_particula.getpz();

    max_px = std::max(max_px, px_dat);
    min_px = std::min(min_px, px_dat);
    max_py = std::max(max_py, py_dat);
    min_py = std::min(min_py, py_dat);
    max_pz = std::max(max_pz, pz_dat);
    min_pz = std::min(min_pz, pz_dat);
  }

  std::vector const max_min{max_px, min_px, max_py, min_py, max_pz, min_pz};
  return max_min;
}
*/
/*
std::vector<double> calculo_tmñ_bloque_malla(std::vector<double> const & max_min,
                                             std::vector<double> const & n_const) {
  double const sx_vect = (max_min[0] - max_min[1]) / n_const[0];
  double const sy_vect = (max_min[2] - max_min[3]) / n_const[1];
  double const sz_vect = (max_min[4] - max_min[5]) / n_const[2];

  return {sx_vect, sy_vect, sz_vect};
}
*/
std::vector<double> longitud_masa(std::ifstream & inputFile) {
  inputFile.seekg(0, std::ios::beg);  // situamos el puntero del archivo en la cabecera
  auto ppm = static_cast<double>(read_binary_value<float>(inputFile));

  // Realizamos los cálculos necesarios
  double const m_dat = constantes::p_const / std::pow(ppm, 3.0);
  double const h_dat = constantes::r_const / ppm;
  return {m_dat, h_dat};
}

int comprobar_fallos_cabecera(std::vector<particula> const & particulas, int n_particulas_int) {
  int const longitud = static_cast<int>(particulas.size());
    int const error = -5;
    if (longitud != n_particulas_int) {
      std::cerr << "Number of particles mismatch. Header:" << n_particulas_int
                << "Found:" << longitud;
      return error;
    }else {
      return 0;
    }
}

// Con esta funcion creamos el array de partículas, y generamos los datos del archivo de entrada
std::pair<int, std::vector<particula>> crear_particulas(std::ifstream & inputFile) {

  inputFile.seekg(4, std::ios::beg);  // situamos el puntero del archivo
  // Obtenemos el número de partículas de la cabecera
  auto n_particulas_int = read_binary_value<int>(inputFile);

  /*TODO:Ese 12 hay que mirarlo*/
  inputFile.seekg(8, std::ios::beg);  // nos aseguramos de empezar después de la cabecera
  std::vector<particula> particulas;  // Creamos el vector con las partículas
  int ident = 0;


  while (!inputFile.eof()) {
    // Leer los datos de una partícula
    int const id_dat = ident;
    auto px_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto py_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto pz_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto hvx         = static_cast<double>(read_binary_value<float>(inputFile));
    auto hvy         = static_cast<double>(read_binary_value<float>(inputFile));
    auto hvz         = static_cast<double>(read_binary_value<float>(inputFile));
    auto vx_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto vy_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto vz_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    particula const n_particula(id_dat, px_dat, py_dat, pz_dat, hvx, hvy, hvz, vy_dat, vx_dat,
                                vz_dat);
    particulas.push_back(n_particula);
    ident++;
  }
  int const error = comprobar_fallos_cabecera(particulas, n_particulas_int);
  if (error < 0 && false) { return std::make_pair(error, std::vector<particula>()); } // este condicional habra q cambiarlo luego
  return std::make_pair(0, particulas);
}

void init_params(std::ifstream & inputFile, int max_iteraciones){
  std::vector<double> vtor =  longitud_masa(inputFile);
  std::pair<int, std::vector<particula>> particulas = crear_particulas(inputFile);
  if (particulas.first == 0){
    using namespace std;
    cout<<"Inizializando malla con m="<<vtor[0]<<" y h="<<vtor[1]<<endl;
    grid malla(vtor, particulas.second);
    using namespace std;
    for(int iteracion = 1; iteracion<=max_iteraciones;iteracion++) {
      cout << "****************************************************" << endl;
      cout << "iniciando  iteracion " << iteracion << endl;

      malla.simular();

      cout<< "finalizada iteracion " << iteracion << endl;
      cout << "----------------------------------------------------" << endl;
    }
  }
}

void particula::colisionLimiteEjeX(double sx, int nx) {
    double const min_value = 0.0000000001;
    int cx = particula::calculoCx(sx, nx);
    if (cx == 0){
        double difLimX = constantes::dp_const - (getpx() - constantes::bmin_const[0]);
        if (difLimX > min_value){
             setax (getax() + (constantes::ps_const * constantes::t_const - constantes::dv_const * getvx()));
        }
    }else if(cx == nx - 1){
        double difLimX = constantes::dp_const - (constantes::bmax_const[0] - getpx());
        if (difLimX > min_value){
            setax(getax() - (constantes::ps_const * difLimX + constantes::dv_const * getvx()));
        }
    }
}

void particula::colisionLimiteEjeY(double sy, int ny) {
    double const min_value = 0.0000000001;
    int cy = particula::calculoCy(sy, ny);
    if (cy == 0){
        double difLimY = constantes::dp_const - (getpy() - constantes::bmin_const[1]);
        if (difLimY > min_value){
            setay(getay() + (constantes::ps_const * constantes::t_const - constantes::dv_const * getvy()));
        }
    }else if(cy == ny - 1){
        double difLimY = constantes::dp_const - (constantes::bmax_const[1] - getpy());
        if (difLimY > min_value){
            setay(getay() - (constantes::ps_const * difLimY + constantes::dv_const * getvy()));
        }
    }
}

void particula::colisionLimiteEjeZ(double sz, int nz) {
    double const min_value = 0.0000000001;
    int cz = particula::calculoCz(sz, nz);
    if (cz == 0){
        double difLimZ = constantes::dp_const - (getpz() - constantes::bmin_const[2]);
        if (difLimZ > min_value){
            setaz(getaz() + (constantes::ps_const * constantes::t_const - constantes::dv_const * getvz()));
        }
    }else if(cz == nz - 1){
        double difLimZ = constantes::dp_const - (constantes::bmax_const[2] - getpz());
        if (difLimZ > min_value){
            setaz(getaz() - (constantes::ps_const * difLimZ + constantes::dv_const * getvz()));
        }
    }
}

int particula::calculoCx(double sx, int nx) {
    double nuevaX = getpx() + gethvx() * constantes::t_const;
    int cx = static_cast<int>(nuevaX / sx);
    if (cx < 0) {
        cx = 0;
    }else if (cx > nx - 1) {
        cx = nx - 1;
    }
    return cx;
}

int particula::calculoCy(double sy, int ny) {
    double nuevaY = getpy() + gethvy() * constantes::t_const;
    int cy = static_cast<int>(nuevaY / sy);
    if (cy < 0) {
        cy = 0;
    }else if (cy > ny - 1) {
        cy = ny - 1;
    }
    return cy;
}

int particula::calculoCz(double sz, int nz) {
    double nuevaZ = getpz() + gethvz() * constantes::t_const;
    int cz = static_cast<int>(nuevaZ / sz);
    if (cz < 0) {
        cz = 0;
    }else if (cz > nz - 1) {
        cz = nz - 1;
    }
    return cz;
}

void particula::actualizarMovimiento() {
    setpx(getpx()+ gethvx() * constantes::t_const + getax() * (constantes::t_const * constantes::t_const));
    setvx(gethvx() + (getax() * constantes::t_const)/2);
    sethvx(gethvx() + getax() * constantes::t_const);

    setpy(getpy()+ gethvy() * constantes::t_const + getay() * (constantes::t_const * constantes::t_const));
    setvy(gethvy() + (getay() * constantes::t_const)/2);
    sethvy(gethvy() + getay() * constantes::t_const);

    setpz(getpz()+ gethvz() * constantes::t_const + getaz() * (constantes::t_const * constantes::t_const));
    setvz(gethvz() + (getaz() * constantes::t_const)/2);
    sethvz(gethvz() + getaz() * constantes::t_const);
}