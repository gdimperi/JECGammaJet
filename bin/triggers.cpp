#include <iostream>
#include <string>
#include <exception>

#include "tinyxml2.h"

#include "triggers.h"

using namespace tinyxml2;

bool Triggers::parse() {
  XMLDocument doc;
  if (doc.LoadFile(mXmlFile.c_str())) {
    doc.PrintError();
    return false;
  }

  const XMLElement* root = doc.FirstChildElement("triggers");
  if (! root)
    return false;

  const XMLElement* runs = root->FirstChildElement("runs");
  for (; runs; runs = runs->NextSiblingElement("runs")) {
    parseRunsElement(runs);
  }

  return true;
}

bool Triggers::parseRunsElement(const XMLElement* runs) {
  Range<unsigned> runRange(runs->UnsignedAttribute("from"), runs->UnsignedAttribute("to"));

  PathVector runPaths;

  const XMLElement* paths = runs->FirstChildElement("path");
  for (; paths; paths = paths->NextSiblingElement("path")) {
    const std::string name = paths->FirstChildElement("name")->GetText();
    const XMLElement* pt = paths->FirstChildElement("pt");

    Range<float> ptRange(pt->FloatAttribute("from"), pt->FloatAttribute("to"));
    runPaths.push_back(std::make_pair(boost::regex(name, boost::regex_constants::icase), ptRange));
  }

  mTriggers[runRange] = runPaths;
  return true;
}

const PathVector& Triggers::getTriggers(unsigned int run) {
  if (mCachedRange && mCachedRange->in(run)) {
    return *mCachedVector;
  }

  for (auto& trigger: mTriggers) {
    const Range<unsigned int>& runRange = trigger.first;

    if (runRange.in(run)) {

      mCachedRange = &runRange;
      mCachedVector = &trigger.second;

      return *mCachedVector;
    }
  }

  throw new std::exception();
}

/*const Regexp& Triggers::getHLTPath(unsigned int run, float pt) {
  const PathVector* paths = NULL;
  if (mCachedRange && mCachedRange->in(run)) {
    paths = mCachedVector;
  } else {
    for (auto& trigger: mTriggers) {
      const Range<unsigned int>& runRange = trigger.first;

      if (runRange.in(run)) {
        mCachedRange = &runRange;
        mCachedVector = &trigger.second;
        paths = mCachedVector;
        break;
      }
    }
  }

  if (! paths) {
    throw new std::exception();
  }

  for (auto& path: *paths) {
    const Range<float>& ptRange = path.second;
    if (ptRange.in(pt)) {
      return path.first;
    }
  }

  throw new std::exception(); // Not found
}*/


void Triggers::print() {
  for (auto& trigger: mTriggers) {
    const Range<unsigned int>& runRange = trigger.first;
    const auto& paths = trigger.second;

    std::cout << "Runs: " << runRange << std::endl;
    for (auto& path: paths) {
      std::cout << path.first << " -> " << path.second << std::endl;
    }
  }
}
