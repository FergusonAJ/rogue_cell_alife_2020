//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once

#include "AbstractGate.h"

class DeterministicGate : public AbstractGate {
public:
  static std::shared_ptr<ParameterLink<std::string>> IO_RangesPL;

  std::vector<std::vector<int>> table;
  DeterministicGate() = delete;
  DeterministicGate(std::shared_ptr<ParametersTable> PT_ = nullptr)
      : AbstractGate(PT_) {
    table = {};
  }
  DeterministicGate(std::pair<std::vector<int>, std::vector<int>> addresses,
                    std::vector<std::vector<int>> _table, int ID_,
                    std::shared_ptr<ParametersTable> PT_ = nullptr);
  virtual ~DeterministicGate() = default;
  virtual void update(std::vector<double> &states,
                      std::vector<double> &nextStates) override;
  virtual std::shared_ptr<AbstractGate>
  makeCopy(std::shared_ptr<ParametersTable> PT_ = nullptr) override;
  // void setupForBits(int* Ins, int nrOfIns, int Out, int logic);
  virtual std::string gateType() override { return "Deterministic"; }
};

