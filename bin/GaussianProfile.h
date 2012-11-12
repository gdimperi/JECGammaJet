#pragma once

#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TDirectory.h>

#include <PhysicsTools/FWLite/interface/TFileService.h>

class GaussianProfile {

  public:
    GaussianProfile(const std::string& name, int nBinsX, const double* binsX, bool doGraph = true):
      m_name(name), m_prefix("pt"), m_autoBinning(true), m_autoBinningLowPercent(0.4), m_autoBinningHighPercent(0.4), m_nXBins(nBinsX), m_XMin(-1), m_XMax(-1), m_dirty(true), m_doGraph(doGraph) {
        m_XBins.assign(binsX, binsX + nBinsX + 1);
      }

    GaussianProfile(const std::string& name, int nBinsX, const double* binsX, int nBinsY, double yMin, double yMax, bool doGraph = true):
      m_name(name), m_prefix("pt"), m_autoBinning(false), m_autoBinningLowPercent(0), m_autoBinningHighPercent(0), m_nXBins(nBinsX), m_XMin(-1), m_XMax(-1),
      m_nYBins(nBinsY), m_YMin(yMin), m_YMax(yMax), m_dirty(true), m_doGraph(doGraph) {
        m_XBins.assign(binsX, binsX + nBinsX + 1);
      }

    GaussianProfile(const std::string& name, int nBinsX, double xMin, double xMax, int nBinsY, double yMin, double yMax, bool doGraph = true):
      m_name(name), m_prefix("pt"), m_autoBinning(false), m_autoBinningLowPercent(0), m_autoBinningHighPercent(0), m_nXBins(nBinsX), m_XMin(xMin), m_XMax(xMax),
      m_nYBins(nBinsY), m_YMin(yMin), m_YMax(yMax), m_dirty(true), m_doGraph(doGraph) {

      }

    void initialize(TFileDirectory& dir) {
      createProfiles(dir);
      mDir = dir.getBareDirectory();
    }

    virtual ~GaussianProfile() {
      /*
      for (TH1* h: m_profiles) {
        delete h;
      }
      */
      write();
    }

    void fill(double x, double y, double weight = 1.0) {
      if (m_profiles.size() == 0) {
        return;
      }

      int bin = findBin(x);
      if (bin < 0)
        return;

      m_profiles[bin]->Fill(y, weight);
      m_dirty = true;
    }

    void drawBin(int bin, Option_t* options) {
      if (bin < 0 || bin >= m_nXBins || m_profiles.size() == 0)
        return;

      m_profiles[bin]->Draw(options);
    }

    void draw(Option_t* options) {
      if (! m_doGraph)
        return;

      if (m_dirty) {
        createGraph();
      }

      m_graph->Draw(options);
    }

    void setAutoBinningPercent(double lowPercent, double highPercent) {
      m_autoBinningLowPercent = lowPercent;
      m_autoBinningHighPercent = highPercent;
    }

    void setPrefix(const std::string prefix) {
      m_prefix = prefix;
    }

    void write() {
      mDir->cd();

      if (m_doGraph && m_dirty) {
        createGraph();
      }

      /*
      for (TH1* h: m_profiles) {
        h->Write();
      }*/

      if (m_doGraph && m_graph.get())
        m_graph->Write();

      /*
      if (f)
        f->cd();
      */
    }

  private:

    void createProfiles(TFileDirectory& dir);
    void createGraph();

    int findBin(double value) const {
      double binWidth = 0;
      if (m_XBins.size() == 0) {
        binWidth = (m_XMax - m_XMin) / (double) m_nXBins;
      }

      for (int i = 0; i < m_nXBins; i++) {

        if (m_XBins.size() > 0) {
          if (value >= m_XBins[i] && value < m_XBins[i + 1])
            return i;
        } else {
          if (value >= (i * binWidth) && value < ((i + 1) * binWidth))
            return i;
        }

      }

      return -1;
    }

    double getBinLowEdge(int bin) const {

      if (bin >= m_nXBins)
        return 0;

      double binWidth = 0;
      if (m_XBins.size() == 0) {
        binWidth = (m_XMax - m_XMin) / (double) m_nXBins;
      }

      if (m_XBins.size() > 0) {
        return m_XBins[bin];
      } else {
        return bin * binWidth;
      }

    }

    double getBinHighEdge(int bin) const {

      if (bin >= m_nXBins)
        return 0;

      double binWidth = 0;
      if (m_XBins.size() == 0) {
        binWidth = (m_XMax - m_XMin) / (double) m_nXBins;
      }

      if (m_XBins.size() > 0) {
        return m_XBins[bin + 1];
      } else {
        return (bin + 1) * binWidth;
      }

    }

    void getTruncatedMeanRMS(TH1* hist, float& mean, float& mean_error, float& rms, float& rms_error);
    void fitProjection(TH1* hist, TF1* gaussian, float nSigma, const std::string& option);

    std::string m_name;
    std::string m_prefix;

    bool m_autoBinning;
    double m_autoBinningLowPercent;
    double m_autoBinningHighPercent;

    int m_nXBins;
    std::vector<double> m_XBins;
    double m_XMin;
    double m_XMax;

    int m_nYBins;
    double m_YMin;
    double m_YMax;

    bool m_dirty;
    std::vector<TH1*> m_profiles;
    std::shared_ptr<TGraphErrors> m_graph;

    bool m_doGraph;

    TDirectory* mDir;
};
