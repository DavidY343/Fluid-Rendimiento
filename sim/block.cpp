//
// Created by david on 11/7/23.
//

#include "block.hpp"

// Función para devolver un vector de partículas que no están dentro del bloque
std::vector<particula> block::devolver_particulas() {
  std::vector<particula> devolver;
  // Iterar a través de las partículas en orden inverso para facilitar la eliminación
  for (int p_dat = static_cast<int>(particulas.size() - 1); 0 <= p_dat; p_dat--) {
    // Verificar si la partícula no está dentro del bloque
    if (!p_bloque(particulas[p_dat])) {
      // Agregar la partícula al vector "devolver" y eliminarla del bloque
      devolver.push_back(particulas[p_dat]);
      this->particulas = eliminar(particulas, p_dat);
    }
  }
  return devolver;
}

// Función para eliminar un elemento en una posición específica de un vector, no deja usar la funcoin erase
std::vector<particula> block::eliminar(std::vector<particula> v, int e) {
  std::vector<particula> eliminado;
  // Iterar sobre el vector original y agregar elementos excepto el elemento en la posición 'e'
  for (unsigned long i = 0; i < v.size(); i++) {
    if (e != static_cast<int>(i)) { eliminado.push_back(v[i]); }
  }
  return eliminado;
}

// Función para añadir una partícula al bloque
void block::anhadir_particulas(particula const & part) {
  particulas.push_back(part);
}

// Función que verifica si una partícula está dentro del bloque
[[nodiscard]] bool block::p_bloque(particula p) const {
  return ((px <= p.getpx()) && (p.getpx() < (px + sx)) && (py <= p.getpy()) &&
          (p.getpy() < (py + sy)) && (pz <= p.getpz()) && (p.getpz() < (pz + sz)));
}
