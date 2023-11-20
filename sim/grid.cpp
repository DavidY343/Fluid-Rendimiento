//
// Created by david on 11/7/23.
//
//
// Created by david on 10/3/23.
//

#include "grid.hpp"

#include "progargs.cpp"
#include "progargs.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>

// Reposiciona las partículas en los bloques según su posición actual
void grid::reposicionar_particulas() {
  using namespace std;
  std::vector<particula> particulas_a_reposicionar;
  for (int indice_bloque = 0; indice_bloque < (nx * ny * nz); indice_bloque++) {
    std::vector<particula> n_part = bloques[indice_bloque].devolver_particulas();
    if (!n_part.empty()) {
      particulas_a_reposicionar.insert(particulas_a_reposicionar.end(), n_part.begin(),
                                       n_part.end());
    }
  }
  for (int par = static_cast<int>(particulas_a_reposicionar.size()) - 1; 0 <= par; par--) {
    recolocar_particula(particulas_a_reposicionar[par]);
  }
}

// Recoloca una partícula en el bloque correspondiente de la malla
void grid::recolocar_particula(particula const & part) {  // esta funcion se puede simplificar
  int indice_i = static_cast<int>((part.getpx() - constantes::bmin_const[0]) / sx);
  if (indice_i < 0) {
    indice_i = 0;
  } else if (indice_i >= nx) {
    indice_i = nx - 1;
  }
  int indice_j = static_cast<int>((part.getpy() - constantes::bmin_const[1]) / sy);
  if (indice_j < 0) {
    indice_j = 0;
  } else if (indice_j >= ny) {
    indice_j = ny - 1;
  }
  int indice_k = static_cast<int>((part.getpz() - constantes::bmin_const[2]) / sz);
  if (indice_k < 0) {
    indice_k = 0;
  } else if (indice_k >= nz) {
    indice_k = nz - 1;
  }
  int const indice_bloque = obtener_indice(indice_i, indice_j, indice_k);
  bloques[indice_bloque].anhadir_particulas(part);
}

// Devuelve los índices de bloques contiguos al bloque especificado
[[nodiscard]] std::vector<int> grid::obtener_contiguos(int n) const {
  std::vector<int> bloques_contiguos;
  std::vector<int> coordenadas = obtener_coordenadas(n);
  for (int i = -1; i < 2; ++i) {
    for (int j = -1; j < 2; ++j) {
      for (int k = -1; k < 2; ++k) {
        if ((i != 0 || j != 0 || k != 0) &&
            (0 <= coordenadas[0] + i && coordenadas[0] + i < nx && 0 <= coordenadas[1] + j &&
             coordenadas[1] + j < ny && 0 <= coordenadas[2] + k && coordenadas[2] + k < nz)) {
          bloques_contiguos.push_back(
              obtener_indice(coordenadas[0] + i, coordenadas[1] + j, coordenadas[2] + k));
        }
      }
    }
  }
  return bloques_contiguos;
}

// Devuelve las coordenadas (índices i, j, k) correspondientes al índice global n
[[nodiscard]] std::vector<int> grid::obtener_coordenadas(int n) const {
  std::vector<int> coordenadas;
  coordenadas.push_back(n / (nz * ny));
  coordenadas.push_back((n % (nz * ny)) / nz);
  coordenadas.push_back((n % (nz * ny)) % nz);
  return coordenadas;
}

// Realiza la simulación, incluyendo el reposicionamiento de partículas, inicialización de densidades,
// cálculo de densidades, transformación de densidades, transferencia de aceleraciones,
// colisiones de partículas y movimiento de partículas
void grid::simular() {
  // 4.3.1
  reposicionar_particulas();
  // 4.3.2
  // inicializar densidades y aceleraciones, a pesar del nombre
  inicializar_densidades();
  calcular_densidades();
  transformar_densidades();
  transferir_aceleraciones();
  // 4.3.3
  colisiones_particulas();
  // 4.3.4 y 4.3.5
  movimiento_particulas();
}

// Inicializa las densidades y aceleraciones de las partículas en cada bloque de la malla
void grid::inicializar_densidades() {
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    for (auto & particula : bloques[indice_bloque].getParticulas()) {
      particula.inicializar_densidad_aceleracion();
    }
  }
}


// Transforma las densidades de las partículas en cada bloque de la malla
void grid::transformar_densidades() {
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    for (auto & particula : bloques[indice_bloque].getParticulas()) {
      particula.transformar_densidad(h, m);
    }
  }
}

/*

  Version en dos bucles
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    for (unsigned long pi = 0; pi < bloques[indice_bloque].getParticulas().size(); pi++) {
      std::vector<int> const bloques_contiguos = obtener_contiguos(indice_bloque);
      for (unsigned long b2 = 0; b2 < bloques_contiguos.size(); b2++) {
        for (unsigned long pj = 0; pj < bloques[b2].getParticulas().size(); pj++) {
          bloques[indice_bloque].getParticulas()[pi].interactuar_aceleracion(bloques[b2].getParticulas()[pj],
  h, m, false);
        }
      }
    }
  }
 */
// Transfiere las aceleraciones entre partiuclas de un mismo bloque y entre particulas de bloque contiguos
// Hemos implementado una fusion de bucle
void grid::transferir_aceleraciones() {
  for (int indice_bloque = 0; indice_bloque < nx * ny * nz; indice_bloque++) {
    std::vector<int> const bloques_contiguos = obtener_contiguos(indice_bloque);
    for (unsigned long pi = 0; pi < bloques[indice_bloque].getParticulas().size(); pi++) {
      for (unsigned long pj = pi + 1; pj < bloques[indice_bloque].getParticulas().size();
           pj++) {  // reducir las iteraciones de este bucle
        bloques[indice_bloque].getParticulas()[pi].interactuar_aceleracion(
            bloques[indice_bloque].getParticulas()[pj], h, m, true);
      }
      // bloques contiguos
      for (int const bloques_contiguo : bloques_contiguos) {
        for (unsigned long pj = 0; pj < bloques[bloques_contiguo].getParticulas().size(); pj++) {
          bloques[indice_bloque].getParticulas()[pi].interactuar_aceleracion(
              bloques[bloques_contiguo].getParticulas()[pj], h, m, false);
        }
      }
    }
  }
}

// Calcula las densidades entre partiuclas de un mismo bloque y entre particulas de bloque contiguos
// Hemos implementado una fusion de bucle
void grid::calcular_densidades() {
  for (int indice_b = 0; indice_b < nx * ny * nz; indice_b++) {
    std::vector<int> const bloques_contiguos = obtener_contiguos(indice_b);
    for (unsigned long pi = 0; pi < bloques[indice_b].getParticulas().size(); pi++) {
      for (unsigned long pj = pi + 1; pj < bloques[indice_b].getParticulas().size(); pj++) {
        bloques[indice_b].getParticulas()[pi].interactuar_densidad(
            bloques[indice_b].getParticulas()[pj], h, true);
      }
      for (int const bloques_contiguo : bloques_contiguos) {
        for (unsigned long pj = 0; pj < bloques[bloques_contiguo].getParticulas().size(); pj++) {
          bloques[indice_b].getParticulas()[pi].interactuar_densidad(
              bloques[bloques_contiguo].getParticulas()[pj], h, false);
        }
      }
    }
  }
}

/*Versión en dos bucles
  for (int b = 0; b < nx * ny * nz; b++) {
    for (unsigned long pi = 0; pi < bloques[b].getParticulas().size(); pi++) {
      std::vector<int> const bloques_contiguos = obtener_contiguos(b);
      for (int const bloques_contiguo : bloques_contiguos) {
        for (unsigned long pj = 0; pj < bloques[bloques_contiguo].getParticulas().size(); pj++) {
          bloques[b].getParticulas()[pi].interactuar_densidad(bloques[bloques_contiguo].getParticulas()[pj],
  h, false);
        }
      }
    }
  }*/

// Realiza colisiones de partículas con los límites de la solo si te encuentras en la esquina de la malla
void grid::colisiones_particulas() {
  std::vector<int> coordenadas;
  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    coordenadas = obtener_coordenadas(i);
    if (coordenadas[0] == 0 || coordenadas[0] == getnx() - 1) {
      bucle_colisiones(i, coordenadas[0] == 0, 0);
    }
    if (coordenadas[1] == 0 || coordenadas[1] == getny() - 1) {
      bucle_colisiones(i, coordenadas[1] == 0, 1);
    }
    if (coordenadas[2] == 0 || coordenadas[2] == getnz() - 1) {
      bucle_colisiones(i, coordenadas[2] == 0, 2);
    }
  }
}

//implementacion del bucle de la funcoin anterior
void grid::bucle_colisiones(int num_bloque, bool lim_inf, int dimension) {
  for (auto & particula : bloques[num_bloque].getParticulas()) {
    switch (dimension) {
      case 0:
        particula.colisionLimiteEjeX(lim_inf);
        break;
      case 1:
        particula.colisionLimiteEjeY(lim_inf);
        break;
      default:
        particula.colisionLimiteEjeZ(lim_inf);
        break;
    }
  }
}

// Actualiza el movimiento de las partículas en la malla
void grid::movimiento_particulas() {
  std::vector<int> coordenadas;
  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    for (auto & particula : bloques[i].getParticulas()) { particula.actualizarMovimiento(); }
    coordenadas = obtener_coordenadas(i);
    if (coordenadas[0] == 0 || coordenadas[0] == getnx() - 1) {
      bucle_limites(i, coordenadas[0] == 0, 0);
    }
    if (coordenadas[1] == 0 || coordenadas[1] == getny() - 1) {
      bucle_limites(i, coordenadas[1] == 0, 1);
    }
    if (coordenadas[2] == 0 || coordenadas[2] == getnz() - 1) {
      bucle_limites(i, coordenadas[2] == 0, 2);
    }
  }
}

void grid::bucle_limites(int num_bloque, bool lim_inf, int dimension) {
  for (auto & particula : bloques[num_bloque].getParticulas()) {
    switch (dimension) {
      case 0:
        particula.limiteRecintox(lim_inf);  // 4.3.5
        break;
      case 1:
        particula.limiteRecintoy(lim_inf);  // 4.3.5
        break;
      default:
        particula.limiteRecintoz(lim_inf);  // 4.3.5
        break;
    }
  }
}

// Imprime los parámetros del sistema y de la malla
void print_params(unsigned long const tamanio, double const ppm, std::vector<double> const & vtor,
                  grid const & malla) {
  using namespace std;
  cout << "Number of particles: " << tamanio << "\n";
  cout << "Particles per meter: " << ppm << "\n";
  cout << "Smoothing length: " << vtor[1] << "\n";
  cout << "Particle mass: " << vtor[0] << "\n";
  cout << "Grid size: " << malla.getnx() << " x " << malla.getny() << " x " << malla.getnz()
       << "\n";
  cout << "Number of blocks: " << malla.getnx() * malla.getny() * malla.getnz() << "\n";
  cout << "Block size: " << malla.getsx() << " x " << malla.getsy() << " x " << malla.getsz()
       << "\n";
}

// Inicia la simulación según los parámetros leidos, es decir, lee el archivo e inicia la simulación
void init_simulation(std::ifstream const & inputFile, int const max_iteraciones,
                     std::ofstream & outputFile) {
  auto ppm           = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
  double const m_dat = constantes::p_const / std::pow(ppm, 3.0);
  double const h_dat = constantes::r_const / ppm;
  std::vector<double> const vtor          = {m_dat, h_dat};
  std::vector<particula> const particulas = crear_particulas(inputFile);
  using namespace std;
  grid malla(vtor, particulas);
  print_params(particulas.size(), ppm, vtor, malla);
  for (int iteracion = 1; iteracion <= max_iteraciones; iteracion++) { malla.simular(); }
  malla.almacenar_resultados(outputFile, ppm, particulas);
}

// Función para convertir los datos de las partículas de doble a simple precisión
std::tuple<float, float, float, float, float, float, float, float, float>
    convertirDatos(particula const & particula) {
  auto px_dat = static_cast<float>(particula.getpx());
  auto py_dat = static_cast<float>(particula.getpy());
  auto pz_dat = static_cast<float>(particula.getpz());
  auto hvx    = static_cast<float>(particula.gethvx());
  auto hvy    = static_cast<float>(particula.gethvy());
  auto hvz    = static_cast<float>(particula.gethvz());
  auto vx_dat = static_cast<float>(particula.getvx());
  auto vy_dat = static_cast<float>(particula.getvy());
  auto vz_dat = static_cast<float>(particula.getvz());

  return std::make_tuple(px_dat, py_dat, pz_dat, hvx, hvy, hvz, vx_dat, vy_dat, vz_dat);
}

// Escribe los datos de una partícula en el archivo de salida
void escribir_datos_particulas(std::ofstream & outputFile, particula const & particula) {
  auto [px_dat, py_dat, pz_dat, hvx, hvy, hvz, vx_dat, vy_dat, vz_dat] = convertirDatos(particula);

  outputFile.write(as_buffer(px_dat), sizeof(px_dat));
  outputFile.write(as_buffer(py_dat), sizeof(py_dat));
  outputFile.write(as_buffer(pz_dat), sizeof(pz_dat));
  outputFile.write(as_buffer(hvx), sizeof(hvx));
  outputFile.write(as_buffer(hvy), sizeof(hvy));
  outputFile.write(as_buffer(hvz), sizeof(hvz));
  outputFile.write(as_buffer(vx_dat), sizeof(vx_dat));
  outputFile.write(as_buffer(vy_dat), sizeof(vy_dat));
  outputFile.write(as_buffer(vz_dat), sizeof(vz_dat));
}

// Función de ordenación personalizada para las partículas
//Solo se ha implementado esto porque no funcian el swap, ademas no de recursion para ordenar
template <typename T>
void mySwap(T & a, T & b) {
  T temp = a;
  a      = b;
  b      = temp;
}

template <typename T>
void mySort(std::vector<T> & arr) {
  for (size_t i = 0; i < arr.size() - 1; ++i) {
    for (size_t j = 0; j < arr.size() - 1 - i; ++j) {
      if (arr[j + 1].getid() < arr[j].getid()) { mySwap(arr[j], arr[j + 1]); }
    }
  }
}

// Almacena los resultados de la simulación en el archivo de salida
void grid::almacenar_resultados(std::ofstream & outputFile, double ppm, std::vector<particula> const & part) {
  outputFile.write(as_buffer(ppm), sizeof(ppm));
  outputFile.write(as_buffer(part.size()), sizeof(part.size()));
  std::vector<particula> particulas;
  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    for (auto & particula : bloques[i].getParticulas()) { particulas.push_back(particula); }
  }
  mySort(particulas);
  // Escribir los datos de las partículas
  for (auto const & particula : particulas) { escribir_datos_particulas(outputFile, particula); }
  outputFile.close();
}

// Metodo para imprimir todas las particulas.
void grid::imprimir_output() {
  std::ofstream outputFile("archivo.out");
  int const precision = 30;  // ajusta según tus necesidades
  outputFile << std::fixed << std::setprecision(precision);
  if (outputFile.is_open()) {
    for (int i = 0; i < getnx() * getny() * getnz(); i++) {
      outputFile << "\n\nBloque " << i << ": " << bloques[i].getParticulas().size()
                 << " particulas: \n";
      for (auto & particula : bloques[i].getParticulas()) {
        outputFile << "Particula " << particula.getid() << ": \n";
        outputFile << "p: " << particula.getpx() << ", " << particula.getpy() << ", "
                   << particula.getpz() << "\n";
        outputFile << "hv: " << particula.gethvx() << ", " << particula.gethvy() << ", "
                   << particula.gethvz() << "\n";
        outputFile << "v: " << particula.getvx() << ", " << particula.getvy() << ", "
                   << particula.getvz() << "\n";
        outputFile << "densidad: " << particula.getp() << "\n";
        outputFile << "a: " << particula.getax() << ", " << particula.getay() << ", "
                   << particula.getaz() << "\n";
      }
    }
    outputFile.close();
    std::cout << "Texto escrito en archivo.out" << '\n';
  } else {
    std::cerr << "No se pudo abrir el archivo.\n";
  }
}

// No me deja usar mergesort pq me dice que es recursivo xD
/*
void merge(std::vector<particula>& arr, size_t l, size_t m, size_t r) {
  size_t n1 = m - l + 1;
  size_t n2 = r - m;
  // Crear vectores temporales
  std::vector<particula> L(arr.begin() + l, arr.begin() + l + n1);
  std::vector<particula> R(arr.begin() + m + 1, arr.begin() + m + 1 + n2);
  // Índices iniciales de los subvectores
  size_t i = 0, j = 0, k = l;
  // Combinar los subvectores de nuevo en arr
  while (i < n1 && j < n2) {
    if (L[i].getid() <= R[j].getid()) {
        arr[k++] = L[i++];
    } else {
        arr[k++] = R[j++];
    }
  }
  // Copiar los elementos restantes de L (si los hay)
  while (i < n1) {
    arr[k++] = L[i++];
  }
  // Copiar los elementos restantes de R (si los hay)
  while (j < n2) {
    arr[k++] = R[j++];
  }
}

// Función principal de merge sort
void mergeSort(std::vector<particula>& arr, size_t l, size_t r) {
  if (l < r) {
    // Encuentra el punto medio
    size_t m = l + (r - l) / 2;

    // Ordena la primera y la segunda mitad
    mergeSort(arr, l, m);
    mergeSort(arr, m + 1, r);

    // Combina las mitades ordenadas
    merge(arr, l, m, r);
  }
}
*/