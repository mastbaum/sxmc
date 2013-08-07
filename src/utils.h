/**
 * \file utils.h
 * \brief Collected utility structures and functions
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLegend.h>

class TNtuple;
class TH1;


/**
 * \struct Range
 * \brief A container for a range of values
 */
template <class T>
struct Range {
  T min;  //!< minimum value
  T max;  //!< maximum value
};


/**
 * Get a value from a TNtuple by event ID and field name.
 *
 * \param nt The source TNtuple
 * \param i The event ID
 * \param field The name of the variable to extract
 * \returns The requested value as a float
 */
float get_ntuple_entry(TNtuple* nt, int i, std::string field);


/**
 * Build a correlation matrix for a TNtuple.
 *
 * Creates a matrix with Pearson product-moment correlation coefficients
 * computed between pairs of variables in a TNtuple. The matrix expressed
 * as a vector of length (entries x entries). Only the upper half is set.
 *
 * \param nt The source TNtuple
 * \returns A correlation matrix as a 1D vector
 */
std::vector<float> get_correlation_matrix(TNtuple* nt);


/**
 * \class SpectralPlot
 * \brief Convenience class for making sensitivity spectral plots
 */
class SpectralPlot {
  public:
    /**
     * Constructor
     * \param _line_width Line width for histograms
     * \param _xmin Minimum of the x range
     * \param _xmax Maximum of the x range
     * \param _ymin Minimum of the y range
     * \param _ymax Maximum of the y range
     * \param _logy Set logscale on the y axis
     * \param _title Histogram title
     * \param _xtitle Title of x axis
     * \param _ytitle Title of y axis
     */
    SpectralPlot(int _line_width=2, float _xmin=1.5, float _xmax=5.0,
                 float _ymin=1e-1, float _ymax=1e3, bool _logy=true,
                 std::string _title="",
                 std::string _xtitle="",
                 std::string _ytitle="Counts/bin");

    /** Copy constructor */
    SpectralPlot(const SpectralPlot& o) {
      this->logy = o.logy;
      this->line_width = o.line_width;
      this->xmin = o.xmin;
      this->xmax = o.xmax;
      this->ymin = o.ymin;
      this->ymax = o.ymax;
      this->title = o.title;
      this->xtitle = o.xtitle;
      this->ytitle = o.ytitle;
      for (size_t i=0; i<o.histograms.size(); i++) {
        this->histograms.push_back(o.histograms[i]);
      }
      this->c = new TCanvas();

      if (o.logy) {
        this->c->SetLogy();
      }
      this->legend = (TLegend*) o.legend->Clone("");
    }

    /** Destructor */
    ~SpectralPlot();

    /**
     * Add a histogram to the plot.
     *
     * No reference to the input histogram is retained.
     *
     * \param _h The histogram
     * \param title Histogram title to show in the legend
     * \param options Options passed to TH1::Draw
     */
    void add(TH1* _h, std::string title, std::string options="");

    /**
     * Save the plot to a file, with type based on extension.
     * 
     * \param filename Name of the output file
     */
    void save(std::string filename);

    /**
     * Create an empty histogram with the same properties as the input.
     *
     * \param h The input histogram
     * \param name Name for the ROOT object
     * \returns An empty histogram, much like the input
     */
    static TH1* make_like(TH1* h, std::string name);

    bool logy;  //!< Log y axis on/off
    int line_width;  //!< Histogram line width
    float xmin;  //!< x-range minimum
    float xmax;  //!< x-range maximum
    float ymin;  //!< y-range minimum
    float ymax;  //!< y-range maximum
    std::string title;  //!< plot title
    std::string xtitle;  //!< x axis title
    std::string ytitle;  //!< y axis title
    TCanvas* c;  //!< Canvas to plot on
    std::vector<TH1*> histograms;  //!< Pointers to plotted histograms
    TLegend* legend;  //!< Plot legend
};

#endif  // __UTILS_H__

