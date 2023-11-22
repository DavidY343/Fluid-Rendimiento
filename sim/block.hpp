#ifndef ARCOS_BLOCK_HPP
#define ARCOS_BLOCK_HPP

#include "progargs.hpp"

#include <iostream>
#include <vector>

class block {
  public:
    // Constructor
    block(std::vector<double> const & position, std::vector<double> const & size)
      : sx(size[0]), sy(size[1]), sz(size[2]), px(position[0]), py(position[1]), pz(position[2]) { }

    // Constructor de copia
    block(block const & other) = default;

    // Operador de asignación de copia
    block & operator=(block const & other) {
      if (this != &other) {  // Check for self-assignment
        sx         = other.sx;
        sy         = other.sy;
        sz         = other.sz;
        px         = other.px;
        py         = other.py;
        pz         = other.pz;
        particulas = other.particulas;
      }
      return *this;
    }

    // Constructor de movimiento
    block(block && other) = default;

    // Destructor
    ~block() = default;

    // Operador de asignación de movimiento
    block & operator=(block && other) noexcept {
      if (this != &other) { }
      return *this;
    }

    std::vector<particula> devolver_particulas();

    static std::vector<particula> eliminar(std::vector<particula> v, int e);

    void anhadir_particulas(particula const & part);

    [[nodiscard]] bool p_bloque(particula p) const;

    [[nodiscard]] std::vector<particula> & getParticulas() { return particulas; }

  private:
    /*Tamaño de bloques*/
    double sx;
    double sy;
    double sz;
    double px;
    double py;
    double pz;
    std::vector<particula> particulas;
};

#endif  // ARCOS_BLOCK_HPP
