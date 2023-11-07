//
// Created by david on 11/7/23.
//

#include "../sim/block.hpp"
#include "../sim/grid.hpp"
#include "../sim/progargs.hpp"

#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>

int main(int argc, char ** argv) {
  std::span const args_view{argv, static_cast<std::size_t>(argc)};
  std::vector<std::string> const argumentos{args_view.begin() + 1, args_view.end()};

  std::ifstream const inputFile(argumentos[1],
                                std::ios::binary);  // Creamos un objeto inputFile para lecutura
  std::ofstream const outputFile(argumentos[2], std::ios::binary);

  int const error = comprobar_params(argc, argumentos, inputFile, outputFile);
  if (error < 0) { return error; }


}