//
// Created by david on 11/7/23.
//

#include "../sim/block.hpp"
#include "../sim/grid.hpp"

#include <fstream>
#include <iostream>
#include <span>
#include <stdexcept>
#include <vector>

int main(int argc, char ** argv) {
  int const nparams = n_params(argc);
  if (nparams < 0) { return nparams; }

  std::span const args_view{argv, static_cast<std::size_t>(argc)};
  std::vector<std::string> const argumentos{args_view.begin() + 1, args_view.end()};

  std::ifstream inputFile(argumentos[1], std::ios::binary);
  std::ofstream outputFile(argumentos[2], std::ios::binary);

  int const error = comprobar_params(argumentos, inputFile, outputFile);
  if (error < 0) { return error; }
  int const max_iteraciones = std::stoi(argumentos[0]);
  init_simulation(inputFile, max_iteraciones, outputFile);
  inputFile.close();
  return (0);
}