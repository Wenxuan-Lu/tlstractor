// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <string>
#include <vector>
#include <memory>
#include <optional>
#include <fstream>
#include <cstring>
#include <zlib.h>
#include <cstdint>
#include <limits>

using namespace Rcpp;

namespace {

// ---------- numerics ----------
static inline double machine_eps() {
  return std::numeric_limits<double>::epsilon(); // ~2.22e-16
}

// Cholesky stability thresholds for genomic reliability
// Two checks: (1) absolute floor catches pathological cases, (2) condition number catches ill-conditioning
// Absolute floor at 1e-6: independent safety regardless of condition number (even n=100 well-conditioned has diagonals ~10)
// Condition number at 1e7: catches relative ill-conditioning (preserves ~8 digits precision)
// Well-conditioned genomic data: typically min_diag ~ 0.1n to n, cond < 100
static constexpr double CHOLESKY_DIAG_FLOOR = 1e-6;         // Absolute minimum (conservative, independent check)
static constexpr double CHOLESKY_COND_LIMIT = 1e7;          // Condition number limit (8-9 digits preserved)
static constexpr double TOL_COLLINEARITY = 1e-14;           // Threshold for detecting x in column space of A
static constexpr double LOGISTIC_BETA_X_FLOOR = 1e-10;      // Threshold for considering beta_x as ~0 for precomp warm start

static inline void set_colpiv_qr_threshold(Eigen::ColPivHouseholderQR<Eigen::MatrixXd>& qr,
                                           const Eigen::Ref<const Eigen::MatrixXd>& X,
                                           double tol_scale = 1.0) {
  // setThreshold internally multiplies the threshold by maxPivot(), so we don't include it here
  const double tol = tol_scale * static_cast<double>(std::max(X.rows(), X.cols()))
    * machine_eps();
  qr.setThreshold(tol);
}

} // namespace

static inline bool ends_with(const std::string& s, const std::string& suf) {
  return s.size() >= suf.size() && s.compare(s.size()-suf.size(), suf.size(), suf) == 0;
}

static inline void trim_cr(std::string& s) {
  if (!s.empty() && s.back() == '\r') s.pop_back();
}

static inline int parse_int_range(const char* p, const char* e) {
  int v = 0;
  while (p < e) { v = v*10 + (*p - '0'); ++p; }
  return v;
}

struct LineReader {
  virtual ~LineReader() = default;
  virtual bool getline(std::string& out) = 0;
};

struct PlainLineReader : LineReader {
  std::ifstream in;
  std::vector<char> buf;
  explicit PlainLineReader(const std::string& path) : in(path.c_str()), buf(1<<20) {
    if (!in) stop("Failed to open input file: " + path);
    in.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
  }
  bool getline(std::string& out) override {
    if (!std::getline(in, out)) return false;
    trim_cr(out);
    return true;
  }
};

struct GzLineReader : LineReader {
  gzFile gz = nullptr;
  std::vector<char> tmp;
  GzLineReader(const std::string& path) : gz(gzopen(path.c_str(), "rb")), tmp(1<<20) {
    if (!gz) stop("Failed to open gz input: " + path);
  }
  ~GzLineReader() override { if (gz) gzclose(gz); }

  bool getline(std::string& out) override {
    out.clear();
    if (!gz) return false;

    while (true) {
      char* res = gzgets(gz, tmp.data(), static_cast<int>(tmp.size()));
      if (!res) {
        if (!out.empty()) {
          trim_cr(out);
          return true;
        }
        return false;
      }
      const size_t len = std::strlen(res);
      out.append(res, len);
      if (!out.empty() && out.back() == '\n') {
        out.pop_back();
        break;
      }
    }
    trim_cr(out);
    return true;
  }
};

struct Writer {
  virtual ~Writer() = default;
  virtual void write(const std::string& s) = 0;
  virtual void flush() = 0;
};

struct PlainWriter : Writer {
  std::ofstream out;
  std::vector<char> buf;
  explicit PlainWriter(const std::string& path)
      : out(path.c_str(), std::ios::out | std::ios::binary),
        buf(1<<20) {
    if (!out) stop("Failed to open output file: " + path);
    out.rdbuf()->pubsetbuf(buf.data(), static_cast<std::streamsize>(buf.size()));
  }
  void write(const std::string& s) override {
    out.write(s.data(), static_cast<std::streamsize>(s.size()));
    if (out.fail()) stop("PlainWriter write failed.");
  }
  void flush() override {
    out.flush();
    if (out.fail()) stop("PlainWriter flush failed.");
  }
};

struct GzWriter : Writer {
  gzFile gz = nullptr;
  explicit GzWriter(const std::string& path)
      : gz(gzopen(path.c_str(), "wb")) { // "wb": Z_DEFAULT_COMPRESSION
    if (!gz) stop("Failed to open gz output: " + path);
    gzbuffer(gz, 1<<20);
  }
  ~GzWriter() override { if (gz) gzclose(gz); }
  void write(const std::string& s) override {
    const char* p = s.data();
    size_t remaining = s.size();
    const size_t max_chunk = static_cast<size_t>(std::numeric_limits<unsigned int>::max());
    while (remaining > 0) {
      const unsigned int chunk = static_cast<unsigned int>(remaining > max_chunk ? max_chunk : remaining);
      const int written = gzwrite(gz, p, chunk);
      if (written != static_cast<int>(chunk)) {
        int errnum = 0;
        const char* msg = gzerror(gz, &errnum);
        stop("gzwrite failed: " + std::string(msg ? msg : "unknown error"));
      }
      p += static_cast<size_t>(chunk);
      remaining -= static_cast<size_t>(chunk);
    }
  }
  void flush() override {
    const int rc = gzflush(gz, Z_SYNC_FLUSH);
    if (rc != Z_OK) {
      int errnum = 0;
      const char* msg = gzerror(gz, &errnum);
      stop("gzflush failed: " + std::string(msg ? msg : "unknown error"));
    }
  }
};

static inline std::unique_ptr<LineReader> open_reader(const std::string& path) {
  if (ends_with(path, ".gz")) return std::unique_ptr<LineReader>(new GzLineReader(path));
  return std::unique_ptr<LineReader>(new PlainLineReader(path));
}
static inline std::unique_ptr<Writer> open_writer(const std::string& path, bool gz) {
  if (gz) return std::unique_ptr<Writer>(new GzWriter(path));
  return std::unique_ptr<Writer>(new PlainWriter(path));
}

// Append small 0..2 values fast
static inline void append_small_0_2(std::string& out, uint8_t v) {
  out.push_back((char)('0' + v));
}

// Find tab positions for first (need_n) fields quickly.
// Returns positions of tabs; tabs[i] is index of i-th tab in the line (0-based).
static inline void find_tabs_firstN(const std::string& line, int need_n_tabs, std::vector<int>& tabs) {
  tabs.clear();
  tabs.reserve(need_n_tabs);
  const char* s = line.c_str();
  int n = (int)line.size();
  for (int i=0; i<n && (int)tabs.size()<need_n_tabs; ++i) {
    if (s[i] == '\t') tabs.push_back(i);
  }
  if ((int)tabs.size() < need_n_tabs) stop("Malformed VCF line or no samples found in header.");
}

// Parse sample "0|1:AN1:AN2" from pointer+len (no allocations)
// Assumes biallelic phased GT with single-digit alleles and no missingness.
static inline void parse_sample_ptr_flare(const char* s, int n, uint8_t& a, uint8_t& b, int& an1, int& an2) {
  a = (uint8_t)(s[0] - '0');
  b = (uint8_t)(s[2] - '0');

  int p = 4;

  an1 = 0;
  while (p < n && s[p] != ':') { an1 = an1*10 + (s[p]-'0'); ++p; }
  ++p;

  an2 = 0;
  while (p < n) { an2 = an2*10 + (s[p]-'0'); ++p; }
}


static inline void parse_gt_ptr(const char* s, int n, uint8_t& a, uint8_t& b) {
  if (n < 3) stop("Malformed GT field.");
  a = (uint8_t)(s[0] - '0');
  b = (uint8_t)(s[2] - '0');
}

static inline void parse_samples_from_txt_tracts_line(const std::string& line, int nsamples, int sample_start, uint8_t* out) {
  const char* s = line.c_str();
  int n = (int)line.size();
  int cur = sample_start;
  for (int i = 0; i < nsamples; ++i) {
    if (cur >= n) stop("Sample count mismatch in text line.");
    out[i] = (uint8_t)(s[cur] - '0');
    cur += 2;
  }
}
static bool read_next_msp_window(LineReader& msp,
                                 std::string& chrom,
                                 int& start,
                                 int& end,
                                 std::vector<int>& calls,
                                 int expected_calls) {
  std::string line;
  while (msp.getline(line)) {
    if (line.empty() || line[0] == '#') continue;

    // Find first 6 tabs to skip to calls
    int tab_count = 0;
    int tpos[6];
    const char* s = line.c_str();
    int n = (int)line.size();
    for (int i=0; i<n && tab_count<6; ++i) {
      if (s[i] == '\t') {
        tpos[tab_count] = i;
        ++tab_count;
      }
    }
    if (tab_count < 6) stop("Malformed MSP line: expected at least 7 columns.");

    chrom.assign(s, s + tpos[0]);
    start = parse_int_range(s + tpos[0] + 1, s + tpos[1]);
    end = parse_int_range(s + tpos[1] + 1, s + tpos[2]);

    calls.clear();
    int cur = tpos[5] + 1;
    while (cur < n) {
      int field_start = cur;
      while (cur < n && s[cur] != '\t') ++cur;
      int field_end = cur;
      if (field_end > field_start) {
        int v = parse_int_range(s + field_start, s + field_end);
        calls.push_back(v);
      }
      ++cur;
    }

    if ((int)calls.size() != expected_calls) {
      stop("MSP calls count mismatch with VCF samples.");
    }
    return true;
  }
  return false;
}

// [[Rcpp::export]]
void extract_tracts_flare_cpp(std::string vcf_path,
                              int num_ancs,
                              std::string out_prefix,
                              CharacterVector formats,
                              int chunk_size,
                              Function gds_append_fn) {
  bool want_txt=false, want_vcf=false, want_gds=false;
  bool gz_txt=false, gz_vcf=false;

  for (int i=0;i<formats.size();++i) {
    std::string f = as<std::string>(formats[i]);
    if (f=="txt") want_txt=true;
    else if (f=="txt.gz") { want_txt=true; gz_txt=true; }
    else if (f=="vcf") want_vcf=true;
    else if (f=="vcf.gz") { want_vcf=true; gz_vcf=true; }
    else if (f=="gds") want_gds=true;
  }

  std::unique_ptr<LineReader> reader = open_reader(vcf_path);

  // Read header lines
  std::string line;
  std::vector<std::string> meta_header;
  std::string col_header;
  std::string txt_header_line;
  int nsamples = 0;

  while (reader->getline(line)) {
    if (line.rfind("##", 0) == 0) { meta_header.push_back(line); continue; }
    if (!line.empty() && line[0] == '#') {
      col_header = line;
      if (!col_header.empty() && col_header.back() == '\t') stop("Malformed VCF header: trailing tab.");
      std::vector<int> tabs;
      find_tabs_firstN(col_header, 9, tabs);
      const int t8 = tabs[8];
      const size_t start = static_cast<size_t>(t8 + 1);
      if (start >= col_header.size()) stop("No samples found in header.");
      int remaining_tabs = 0;
      for (size_t i = start; i < col_header.size(); ++i) {
        if (col_header[i] == '\t') ++remaining_tabs;
      }
      nsamples = remaining_tabs + 1;
      if (nsamples <= 1) stop("Error: at least two samples need to be present.");
      txt_header_line = "CHROM\tPOS\tID\tREF\tALT\t" + col_header.substr(static_cast<size_t>(t8 + 1)) + "\n";
      break;
    }
  }
  if (col_header.empty()) stop("Failed to find VCF header (#CHROM line).");
 
  // Writers
  std::vector< std::unique_ptr<Writer> > dos_w(num_ancs), hap_w(num_ancs), vcf_w(num_ancs);

  if (want_txt) {
    for (int k=0;k<num_ancs;++k) {
      dos_w[k] = open_writer(out_prefix + ".anc" + std::to_string(k) + ".dosage.txt" + (gz_txt ? ".gz" : ""), gz_txt);
      hap_w[k] = open_writer(out_prefix + ".anc" + std::to_string(k) + ".hapcount.txt" + (gz_txt ? ".gz" : ""), gz_txt);
      dos_w[k]->write(txt_header_line);
      hap_w[k]->write(txt_header_line);
    }
  }

  if (want_vcf) {
    const std::string newline("\n");
    for (int k=0;k<num_ancs;++k) {
      vcf_w[k] = open_writer(out_prefix + ".anc" + std::to_string(k) + ".vcf" + (gz_vcf ? ".gz" : ""), gz_vcf);
      for (const auto &mh : meta_header) {
        vcf_w[k]->write(mh);
        vcf_w[k]->write(newline);
      }
      vcf_w[k]->write(col_header);
      vcf_w[k]->write(newline);
    }
  }

  std::vector<std::string>().swap(meta_header);
  std::string().swap(col_header);
  std::string().swap(txt_header_line);

  // Working arrays (flat: [anc * nsamples + i])
  const size_t total = (size_t)num_ancs * (size_t)nsamples;
  std::vector<uint8_t> dos_cur(total), hap_cur(total);
  std::vector<uint8_t> a(nsamples), b(nsamples);
  std::vector<int> an1(nsamples), an2(nsamples);

  // Output buffers for TXT/VCF
  std::vector<std::string> buf_dos(num_ancs), buf_hap(num_ancs), buf_vcf(num_ancs);
  const size_t meta_est = (size_t)chunk_size * 64; // rough per-line metadata estimate
  for (int k=0;k<num_ancs;++k) {
    if (want_txt) {
      buf_dos[k].reserve(meta_est + (size_t)chunk_size * (size_t)nsamples * 2);
      buf_hap[k].reserve(meta_est + (size_t)chunk_size * (size_t)nsamples * 2);
    }
    if (want_vcf) {
      buf_vcf[k].reserve(meta_est + (size_t)chunk_size * (size_t)nsamples * 4);
    }
  }

  // GDS chunk buffers (raw int8 matrices per ancestry)
  std::vector< RawVector > gds_dos(num_ancs), gds_hap(num_ancs);
  CharacterVector r_chrom, r_id, r_ref, r_alt;
  IntegerVector r_pos;
  if (want_gds) {
    for (int k=0;k<num_ancs;++k) {
      gds_dos[k] = RawVector((R_xlen_t)nsamples * (R_xlen_t)chunk_size);
      gds_hap[k] = RawVector((R_xlen_t)nsamples * (R_xlen_t)chunk_size);
      gds_dos[k].attr("dim") = IntegerVector::create(nsamples, chunk_size);
      gds_hap[k].attr("dim") = IntegerVector::create(nsamples, chunk_size);
    }
    r_chrom = CharacterVector(chunk_size);
    r_id    = CharacterVector(chunk_size);
    r_ref   = CharacterVector(chunk_size);
    r_alt   = CharacterVector(chunk_size);
    r_pos   = IntegerVector(chunk_size);
  }

  // Tabs buffer for extracting fixed fields and sample boundaries
  std::vector<int> tabs9;
  int since_flush = 0;

  // Process variant lines
  while (reader->getline(line)) {
    if (line.empty() || line[0] == '#') continue;

    // Full reset per variant (memset)
    std::memset(dos_cur.data(), 0, total * sizeof(uint8_t));
    std::memset(hap_cur.data(), 0, total * sizeof(uint8_t));

    // Locate first 9 tabs (0..8) to define fixed columns and start of sample fields
    find_tabs_firstN(line, 9, tabs9);
    const char* s = line.c_str();
    const int n = (int)line.size();

    int t0 = tabs9[0], t1=tabs9[1], t2=tabs9[2], t3=tabs9[3], t4=tabs9[4], t5=tabs9[5], t6=tabs9[6], t7=tabs9[7], t8=tabs9[8];

    const char* chrom_p = s;
    const char* pos_p = s + t0 + 1;
    const char* pos_e = s + t1;
    int posv = parse_int_range(pos_p, pos_e);
    const char* id_p = s + t1 + 1;
    const char* ref_p = s + t2 + 1;
    const char* alt_p = s + t3 + 1;
    
    // Parse samples without allocating substrings:
    // sample fields start at t8+1
    int col_start = t8 + 1;
    int sample_i = 0;
    int cur = col_start;

    while (cur <= n && sample_i < nsamples) {
      int start = cur;
      while (cur < n && s[cur] != '\t') ++cur;
      int end = cur;
      // parse sample field [start, end)
      uint8_t ai, bi;
      int an1i, an2i;
      parse_sample_ptr_flare(s + start, end - start, ai, bi, an1i, an2i);

      a[sample_i] = ai;
      b[sample_i] = bi;
      an1[sample_i] = an1i;
      an2[sample_i] = an2i;

      if (an1i < 0 || an1i >= num_ancs || an2i < 0 || an2i >= num_ancs) stop("Ancestry code out of range.");

      // update dosage/hapcount at (an1, sample_i) and (an2, sample_i)
      size_t idx1 = (size_t)an1i * (size_t)nsamples + (size_t)sample_i;
      size_t idx2 = (size_t)an2i * (size_t)nsamples + (size_t)sample_i;

      hap_cur[idx1] = (uint8_t)(hap_cur[idx1] + 1);
      hap_cur[idx2] = (uint8_t)(hap_cur[idx2] + 1);
      dos_cur[idx1] = (uint8_t)(dos_cur[idx1] + ai);
      dos_cur[idx2] = (uint8_t)(dos_cur[idx2] + bi);

      ++sample_i;
      ++cur; // skip tab
    }

    if (sample_i != nsamples) stop("Sample count mismatch when parsing variant line.");

    // TXT output
    if (want_txt) {
      for (int k=0;k<num_ancs;++k) {
        size_t base = (size_t)k * (size_t)nsamples;

        // dosage line
        std::string& od = buf_dos[k];
        od.append(s, s + t4);

        for (int i=0;i<nsamples;++i) {
          od.push_back('\t');
          append_small_0_2(od, dos_cur[base + (size_t)i]);
        }
        od.push_back('\n');

        // hapcount line
        std::string& oh = buf_hap[k];
        oh.append(s, s + t4);

        for (int i=0;i<nsamples;++i) {
          oh.push_back('\t');
          append_small_0_2(oh, hap_cur[base + (size_t)i]);
        }
        oh.push_back('\n');
      }
    }

    // VCF output
    if (want_vcf) {
      std::string prefix;
      prefix.reserve((size_t)t7 + 16);
      prefix.append(s, s + t7);
      prefix.append("\tGT");

      for (int k=0;k<num_ancs;++k) {
        std::string& ov = buf_vcf[k];
        ov.append(prefix);
        for (int i=0;i<nsamples;++i) {
          ov.push_back('\t');
          char h1 = (an1[i] == k) ? (char)('0' + a[i]) : '.';
          char h2 = (an2[i] == k) ? (char)('0' + b[i]) : '.';
          ov.push_back(h1);
          ov.push_back('|');
          ov.push_back(h2);
        }
        ov.push_back('\n');
      }
    }

    // GDS chunk fill (use memcpy row slices per ancestry)
    if (want_gds) {
      SET_STRING_ELT(r_chrom, since_flush, Rf_mkCharLenCE(chrom_p, t0, CE_UTF8));
      SET_STRING_ELT(r_id, since_flush, Rf_mkCharLenCE(id_p, t2 - t1 - 1, CE_UTF8));
      SET_STRING_ELT(r_ref, since_flush, Rf_mkCharLenCE(ref_p, t3 - t2 - 1, CE_UTF8)); 
      SET_STRING_ELT(r_alt, since_flush, Rf_mkCharLenCE(alt_p, t4 - t3 - 1, CE_UTF8));
      r_pos[since_flush] = posv;

      R_xlen_t col_off = (R_xlen_t)since_flush * (R_xlen_t)nsamples;
      for (int k=0;k<num_ancs;++k) {
        const uint8_t* src_d = dos_cur.data() + (size_t)k * (size_t)nsamples;
        const uint8_t* src_h = hap_cur.data() + (size_t)k * (size_t)nsamples;
        std::memcpy(&gds_dos[k][col_off], src_d, (size_t)nsamples * sizeof(uint8_t));
        std::memcpy(&gds_hap[k][col_off], src_h, (size_t)nsamples * sizeof(uint8_t));
      }
    }

    // Flush TXT/VCF buffers every chunk_size variants (reduces syscalls)
    ++since_flush;
    if (since_flush >= chunk_size) {
      if (want_txt) {
        for (int k=0;k<num_ancs;++k) {
          dos_w[k]->write(buf_dos[k]); buf_dos[k].clear();
          hap_w[k]->write(buf_hap[k]); buf_hap[k].clear();
        }
      }
      if (want_vcf) {
        for (int k=0;k<num_ancs;++k) {
          vcf_w[k]->write(buf_vcf[k]); buf_vcf[k].clear();
        }
      }
      if (want_gds) {
        List dos_list(num_ancs), hap_list(num_ancs);
        for (int k=0;k<num_ancs;++k) {
          dos_list[k] = gds_dos[k];
          hap_list[k] = gds_hap[k];
        }
        gds_append_fn(r_chrom, r_pos, r_id, r_ref, r_alt, dos_list, hap_list);
      }
      since_flush = 0;
    }
  }

  // Final flush
  if (want_txt) {
    for (int k=0;k<num_ancs;++k) {
      if (!buf_dos[k].empty()) dos_w[k]->write(buf_dos[k]);
      if (!buf_hap[k].empty()) hap_w[k]->write(buf_hap[k]);
      dos_w[k]->flush();
      hap_w[k]->flush();
    }
  }
  if (want_vcf) {
    for (int k=0;k<num_ancs;++k) {
      if (!buf_vcf[k].empty()) vcf_w[k]->write(buf_vcf[k]);
      vcf_w[k]->flush();
    }
  }
  if (want_gds && since_flush > 0) {
    List dos_list(num_ancs), hap_list(num_ancs);
    for (int k=0;k<num_ancs;++k) {
      const R_xlen_t ncopy = (R_xlen_t)nsamples * (R_xlen_t)since_flush;
      RawVector dos_part(ncopy);
      RawVector hap_part(ncopy);
      std::memcpy(&dos_part[0], &gds_dos[k][0], (size_t)ncopy * sizeof(uint8_t));
      std::memcpy(&hap_part[0], &gds_hap[k][0], (size_t)ncopy * sizeof(uint8_t));
      dos_part.attr("dim") = IntegerVector::create(nsamples, since_flush);
      hap_part.attr("dim") = IntegerVector::create(nsamples, since_flush);
      dos_list[k] = dos_part;
      hap_list[k] = hap_part;
    }
    CharacterVector rchrom = r_chrom[Range(0, since_flush - 1)];
    CharacterVector rid    = r_id[Range(0, since_flush - 1)];
    CharacterVector rref   = r_ref[Range(0, since_flush - 1)];
    CharacterVector ralt   = r_alt[Range(0, since_flush - 1)];
    IntegerVector rpos     = r_pos[Range(0, since_flush - 1)];
    gds_append_fn(rchrom, rpos, rid, rref, ralt, dos_list, hap_list);
  }

  return;
}

// [[Rcpp::export]]
void extract_tracts_cpp(std::string vcf_path,
                        std::string msp_path,
                        int num_ancs,
                        std::string out_prefix,
                        CharacterVector formats,
                        int chunk_size,
                        Function gds_append_fn) {
  bool want_txt=false, want_vcf=false, want_gds=false;
  bool gz_txt=false, gz_vcf=false;

  for (int i=0;i<formats.size();++i) {
    std::string f = as<std::string>(formats[i]);
    if (f=="txt") want_txt=true;
    else if (f=="txt.gz") { want_txt=true; gz_txt=true; }
    else if (f=="vcf") want_vcf=true;
    else if (f=="vcf.gz") { want_vcf=true; gz_vcf=true; }
    else if (f=="gds") want_gds=true;
  }

  std::unique_ptr<LineReader> vcf = open_reader(vcf_path);
  std::unique_ptr<LineReader> msp = open_reader(msp_path);

  // Read VCF header
  std::string line;
  std::vector<std::string> meta_header;
  std::string col_header;
  std::string txt_header_line;
  int nsamples = 0;

  while (vcf->getline(line)) {
    if (line.rfind("##", 0) == 0) {
      if (want_vcf) meta_header.push_back(line);
      continue;
    }
    if (!line.empty() && line[0] == '#') {
      col_header = line;
      if (!col_header.empty() && col_header.back() == '\t') stop("Malformed VCF header: trailing tab.");
      std::vector<int> tabs;
      find_tabs_firstN(col_header, 9, tabs);
      const int t8 = tabs[8];
      const size_t start = static_cast<size_t>(t8 + 1);
      if (start >= col_header.size()) stop("No samples found in header.");
      int remaining_tabs = 0;
      for (size_t i = start; i < col_header.size(); ++i) {
        if (col_header[i] == '\t') ++remaining_tabs;
      }
      nsamples = remaining_tabs + 1;
      if (nsamples <= 1) stop("Error: at least two samples need to be present.");
      txt_header_line = "CHROM\tPOS\tID\tREF\tALT\t" + col_header.substr(static_cast<size_t>(t8 + 1)) + "\n";
      break;
    }
  }
  if (col_header.empty()) stop("Failed to find VCF header (#CHROM line).");

  // Writers
  std::vector< std::unique_ptr<Writer> > dos_w(num_ancs), hap_w(num_ancs), vcf_w(num_ancs);
  if (want_txt) {
    for (int k=0;k<num_ancs;++k) {
      dos_w[k] = open_writer(out_prefix + ".anc" + std::to_string(k) + ".dosage.txt" + (gz_txt ? ".gz" : ""), gz_txt);
      hap_w[k] = open_writer(out_prefix + ".anc" + std::to_string(k) + ".hapcount.txt" + (gz_txt ? ".gz" : ""), gz_txt);
      dos_w[k]->write(txt_header_line);
      hap_w[k]->write(txt_header_line);
    }
  }

  if (want_vcf) {
    const std::string newline("\n");
    for (int k=0;k<num_ancs;++k) {
      vcf_w[k] = open_writer(out_prefix + ".anc" + std::to_string(k) + ".vcf" + (gz_vcf ? ".gz" : ""), gz_vcf);
      for (const auto &mh : meta_header) {
        vcf_w[k]->write(mh);
        vcf_w[k]->write(newline);
      }
      vcf_w[k]->write(col_header);
      vcf_w[k]->write(newline);
    }
  }

  std::vector<std::string>().swap(meta_header);
  std::string().swap(col_header);
  std::string().swap(txt_header_line);

  const size_t total = (size_t)num_ancs * (size_t)nsamples;
  std::vector<uint8_t> dos_cur(total), hap_cur(total);
  std::vector<uint8_t> a(nsamples), b(nsamples);
  std::vector<uint8_t> an1(nsamples), an2(nsamples);

  std::vector<std::string> buf_dos(num_ancs), buf_hap(num_ancs), buf_vcf(num_ancs);
  const size_t meta_est = (size_t)chunk_size * 64;
  for (int k=0;k<num_ancs;++k) {
    if (want_txt) {
      buf_dos[k].reserve(meta_est + (size_t)chunk_size * (size_t)nsamples * 2);
      buf_hap[k].reserve(meta_est + (size_t)chunk_size * (size_t)nsamples * 2);
    }
    if (want_vcf) {
      buf_vcf[k].reserve(meta_est + (size_t)chunk_size * (size_t)nsamples * 4);
    }
  }

  // GDS chunk buffers
  std::vector< RawVector > gds_dos(num_ancs), gds_hap(num_ancs);
  CharacterVector r_chrom, r_id, r_ref, r_alt;
  IntegerVector r_pos;
  if (want_gds) {
    for (int k=0;k<num_ancs;++k) {
      gds_dos[k] = RawVector((R_xlen_t)nsamples * (R_xlen_t)chunk_size);
      gds_hap[k] = RawVector((R_xlen_t)nsamples * (R_xlen_t)chunk_size);
      gds_dos[k].attr("dim") = IntegerVector::create(nsamples, chunk_size);
      gds_hap[k].attr("dim") = IntegerVector::create(nsamples, chunk_size);
    }
    r_chrom = CharacterVector(chunk_size);
    r_id    = CharacterVector(chunk_size);
    r_ref   = CharacterVector(chunk_size);
    r_alt   = CharacterVector(chunk_size);
    r_pos   = IntegerVector(chunk_size);
  }

  std::vector<int> tabs9;
  int since_flush = 0;

  std::string win_chrom;
  int win_start = 0, win_end = 0;
  bool have_window = false;
  std::vector<int> calls;
  calls.reserve((size_t)nsamples * 2);

  while (vcf->getline(line)) {
    if (line.empty() || line[0] == '#') continue;

    find_tabs_firstN(line, 9, tabs9);
    const char* s = line.c_str();
    const int n = (int)line.size();

    int t0 = tabs9[0], t1=tabs9[1], t2=tabs9[2], t3=tabs9[3], t4=tabs9[4], t5=tabs9[5], t6=tabs9[6], t7=tabs9[7], t8=tabs9[8];

    std::string vcf_chrom(s, s + t0);
    const char* pos_p = s + t0 + 1;
    const char* pos_e = s + t1;
    int posv = parse_int_range(pos_p, pos_e);
    const char* id_p = s + t1 + 1;
    const char* ref_p = s + t2 + 1;
    const char* alt_p = s + t3 + 1;

    bool skip_line = false;
    while (!(have_window && vcf_chrom == win_chrom && win_start <= posv && posv < win_end)) {
      if (have_window && vcf_chrom == win_chrom && win_start > posv) {
        skip_line = true;
        break;
      }
      have_window = read_next_msp_window(*msp, win_chrom, win_start, win_end, calls, nsamples * 2);
      if (!have_window) break;
      if (vcf_chrom == win_chrom && win_start > posv) {
        skip_line = true;
        break;
      }
      if (vcf_chrom != win_chrom) {
        stop("VCF chromosome does not match MSP chromosome.");
      }
    }
    if (!have_window) break;
    if (skip_line) continue;

    std::memset(dos_cur.data(), 0, total * sizeof(uint8_t));
    std::memset(hap_cur.data(), 0, total * sizeof(uint8_t));

    int col_start = t8 + 1;
    int sample_i = 0;
    int cur = col_start;
    while (cur <= n && sample_i < nsamples) {
      int start = cur;
      while (cur < n && s[cur] != '\t') ++cur;
      int end = cur;

      uint8_t ai, bi;
      parse_gt_ptr(s + start, end - start, ai, bi);
      a[sample_i] = ai;
      b[sample_i] = bi;

      int anc_a = calls[2*sample_i];
      int anc_b = calls[2*sample_i + 1];
      if (anc_a < 0 || anc_a >= num_ancs || anc_b < 0 || anc_b >= num_ancs) stop("Ancestry code out of range.");

      an1[sample_i] = (uint8_t)anc_a;
      an2[sample_i] = (uint8_t)anc_b;

      size_t idx1 = (size_t)anc_a * (size_t)nsamples + (size_t)sample_i;
      size_t idx2 = (size_t)anc_b * (size_t)nsamples + (size_t)sample_i;

      hap_cur[idx1] = (uint8_t)(hap_cur[idx1] + 1);
      hap_cur[idx2] = (uint8_t)(hap_cur[idx2] + 1);
      dos_cur[idx1] = (uint8_t)(dos_cur[idx1] + ai);
      dos_cur[idx2] = (uint8_t)(dos_cur[idx2] + bi);

      ++sample_i;
      ++cur;
    }
    if (sample_i != nsamples) stop("Sample count mismatch when parsing VCF line.");

    if (want_txt) {
      for (int k=0;k<num_ancs;++k) {
        size_t base = (size_t)k * (size_t)nsamples;
        std::string& od = buf_dos[k];
        od.append(s, s + t4);
        for (int i=0;i<nsamples;++i) {
          od.push_back('\t');
          append_small_0_2(od, dos_cur[base + (size_t)i]);
        }
        od.push_back('\n');

        std::string& oh = buf_hap[k];
        oh.append(s, s + t4);
        for (int i=0;i<nsamples;++i) {
          oh.push_back('\t');
          append_small_0_2(oh, hap_cur[base + (size_t)i]);
        }
        oh.push_back('\n');
      }
    }

    if (want_vcf) {
      std::string prefix;
      prefix.reserve((size_t)t7 + 16);
      prefix.append(s, s + t7);
      prefix.append("\tGT");

      for (int k=0;k<num_ancs;++k) {
        std::string& ov = buf_vcf[k];
        ov.append(prefix);
        for (int i=0;i<nsamples;++i) {
          ov.push_back('\t');
          char h1 = (an1[i] == k) ? (char)('0' + a[i]) : '.';
          char h2 = (an2[i] == k) ? (char)('0' + b[i]) : '.';
          ov.push_back(h1);
          ov.push_back('|');
          ov.push_back(h2);
        }
        ov.push_back('\n');
      }
    }

    if (want_gds) {
      SET_STRING_ELT(r_chrom, since_flush, Rf_mkCharLenCE(s, t0, CE_UTF8));
      SET_STRING_ELT(r_id, since_flush, Rf_mkCharLenCE(id_p, t2 - t1 - 1, CE_UTF8));
      SET_STRING_ELT(r_ref, since_flush, Rf_mkCharLenCE(ref_p, t3 - t2 - 1, CE_UTF8));
      SET_STRING_ELT(r_alt, since_flush, Rf_mkCharLenCE(alt_p, t4 - t3 - 1, CE_UTF8));
      r_pos[since_flush] = posv;

      R_xlen_t col_off = (R_xlen_t)since_flush * (R_xlen_t)nsamples;
      for (int k=0;k<num_ancs;++k) {
        const uint8_t* src_d = dos_cur.data() + (size_t)k * (size_t)nsamples;
        const uint8_t* src_h = hap_cur.data() + (size_t)k * (size_t)nsamples;
        std::memcpy(&gds_dos[k][col_off], src_d, (size_t)nsamples * sizeof(uint8_t));
        std::memcpy(&gds_hap[k][col_off], src_h, (size_t)nsamples * sizeof(uint8_t));
      }
    }

    ++since_flush;
    if (since_flush >= chunk_size) {
      if (want_txt) {
        for (int k=0;k<num_ancs;++k) {
          dos_w[k]->write(buf_dos[k]); buf_dos[k].clear();
          hap_w[k]->write(buf_hap[k]); buf_hap[k].clear();
        }
      }
      if (want_vcf) {
        for (int k=0;k<num_ancs;++k) {
          vcf_w[k]->write(buf_vcf[k]); buf_vcf[k].clear();
        }
      }
      if (want_gds) {
        List dos_list(num_ancs), hap_list(num_ancs);
        for (int k=0;k<num_ancs;++k) {
          dos_list[k] = gds_dos[k];
          hap_list[k] = gds_hap[k];
        }
        gds_append_fn(r_chrom, r_pos, r_id, r_ref, r_alt, dos_list, hap_list);
      }
      since_flush = 0;
    }
  }

  if (want_txt) {
    for (int k=0;k<num_ancs;++k) {
      if (!buf_dos[k].empty()) dos_w[k]->write(buf_dos[k]);
      if (!buf_hap[k].empty()) hap_w[k]->write(buf_hap[k]);
      dos_w[k]->flush();
      hap_w[k]->flush();
    }
  }
  if (want_vcf) {
    for (int k=0;k<num_ancs;++k) {
      if (!buf_vcf[k].empty()) vcf_w[k]->write(buf_vcf[k]);
      vcf_w[k]->flush();
    }
  }
  if (want_gds && since_flush > 0) {
    List dos_list(num_ancs), hap_list(num_ancs);
    for (int k=0;k<num_ancs;++k) {
      const R_xlen_t ncopy = (R_xlen_t)nsamples * (R_xlen_t)since_flush;
      RawVector dos_part(ncopy);
      RawVector hap_part(ncopy);
      std::memcpy(&dos_part[0], &gds_dos[k][0], (size_t)ncopy * sizeof(uint8_t));
      std::memcpy(&hap_part[0], &gds_hap[k][0], (size_t)ncopy * sizeof(uint8_t));
      dos_part.attr("dim") = IntegerVector::create(nsamples, since_flush);
      hap_part.attr("dim") = IntegerVector::create(nsamples, since_flush);
      dos_list[k] = dos_part;
      hap_list[k] = hap_part;
    }
    CharacterVector rchrom = r_chrom[Range(0, since_flush - 1)];
    CharacterVector rid    = r_id[Range(0, since_flush - 1)];
    CharacterVector rref   = r_ref[Range(0, since_flush - 1)];
    CharacterVector ralt   = r_alt[Range(0, since_flush - 1)];
    IntegerVector rpos     = r_pos[Range(0, since_flush - 1)];
    gds_append_fn(rchrom, rpos, rid, rref, ralt, dos_list, hap_list);
  }

  return;
}

// [[Rcpp::export]]
void convert_tracts_txt_to_gds_cpp(CharacterVector dosage_files,
                                  CharacterVector hap_files,
                                  int nsamples,
                                  int chunk_size,
                                  Function gds_append_fn) {
  const int num_ancs = dosage_files.size();
  
  std::vector<std::unique_ptr<LineReader>> dos_r(num_ancs), hap_r(num_ancs);
  for (int k=0;k<num_ancs;++k) {
    dos_r[k] = open_reader(as<std::string>(dosage_files[k]));
    hap_r[k] = open_reader(as<std::string>(hap_files[k]));
  }

  std::string line;
  for (int k=0;k<num_ancs;++k) {
    if (!dos_r[k]->getline(line)) stop("Failed to read dosage header line.");
    if (!hap_r[k]->getline(line)) stop("Failed to read hapcount header line.");
  }

  std::vector<uint8_t> dos_cur((size_t)num_ancs * (size_t)nsamples);
  std::vector<uint8_t> hap_cur((size_t)num_ancs * (size_t)nsamples);

  std::vector< RawVector > gds_dos(num_ancs), gds_hap(num_ancs);
  CharacterVector r_chrom(chunk_size), r_id(chunk_size), r_ref(chunk_size), r_alt(chunk_size);
  IntegerVector r_pos(chunk_size);
  for (int k=0;k<num_ancs;++k) {
    gds_dos[k] = RawVector((R_xlen_t)nsamples * (R_xlen_t)chunk_size);
    gds_hap[k] = RawVector((R_xlen_t)nsamples * (R_xlen_t)chunk_size);
    gds_dos[k].attr("dim") = IntegerVector::create(nsamples, chunk_size);
    gds_hap[k].attr("dim") = IntegerVector::create(nsamples, chunk_size);
  }

  std::vector<std::string> dos_line(num_ancs), hap_line(num_ancs);
  std::vector<int> tabs5;
  int since_flush = 0;

  while (true) {
    if (!dos_r[0]->getline(dos_line[0])) break;
    for (int k=1;k<num_ancs;++k) {
      if (!dos_r[k]->getline(dos_line[k])) stop("Dosage file line count mismatch.");
    }
    for (int k=0;k<num_ancs;++k) {
      if (!hap_r[k]->getline(hap_line[k])) stop("Hapcount file line count mismatch.");
    }

    find_tabs_firstN(dos_line[0], 5, tabs5);
    const char* s = dos_line[0].c_str();
    int t0 = tabs5[0], t1 = tabs5[1], t2 = tabs5[2], t3 = tabs5[3], t4 = tabs5[4];
    const char* chrom_p = s;
    const char* pos_p = s + t0 + 1;
    const char* pos_e = s + t1;
    int posv = parse_int_range(pos_p, pos_e);
    const char* id_p = s + t1 + 1;
    const char* ref_p = s + t2 + 1;
    const char* alt_p = s + t3 + 1;

    SET_STRING_ELT(r_chrom, since_flush, Rf_mkCharLenCE(chrom_p, t0, CE_UTF8));
    SET_STRING_ELT(r_id, since_flush, Rf_mkCharLenCE(id_p, t2 - t1 - 1, CE_UTF8));
    SET_STRING_ELT(r_ref, since_flush, Rf_mkCharLenCE(ref_p, t3 - t2 - 1, CE_UTF8));
    SET_STRING_ELT(r_alt, since_flush, Rf_mkCharLenCE(alt_p, t4 - t3 - 1, CE_UTF8));
    r_pos[since_flush] = posv;

    for (int k=0;k<num_ancs;++k) {
      const size_t k_off = (size_t)k * (size_t)nsamples;
      parse_samples_from_txt_tracts_line(dos_line[k], nsamples, t4 + 1, dos_cur.data() + k_off);
      parse_samples_from_txt_tracts_line(hap_line[k], nsamples, t4 + 1, hap_cur.data() + k_off);
    }

    R_xlen_t col_off = (R_xlen_t)since_flush * (R_xlen_t)nsamples;
    const size_t byte_count = (size_t)nsamples * sizeof(uint8_t);
    for (int k=0;k<num_ancs;++k) {
      const size_t k_off = (size_t)k * (size_t)nsamples;
      const uint8_t* src_d = dos_cur.data() + k_off;
      const uint8_t* src_h = hap_cur.data() + k_off;
      std::memcpy(&gds_dos[k][col_off], src_d, byte_count);
      std::memcpy(&gds_hap[k][col_off], src_h, byte_count);
    }

    ++since_flush;
    if (since_flush >= chunk_size) {
      List dos_list(num_ancs), hap_list(num_ancs);
      for (int k=0;k<num_ancs;++k) {
        dos_list[k] = gds_dos[k];
        hap_list[k] = gds_hap[k];
      }
      gds_append_fn(r_chrom, r_pos, r_id, r_ref, r_alt, dos_list, hap_list);
      since_flush = 0;
    }
  }

  if (since_flush > 0) {
    List dos_list(num_ancs), hap_list(num_ancs);
    for (int k=0;k<num_ancs;++k) {
      const R_xlen_t ncopy = (R_xlen_t)nsamples * (R_xlen_t)since_flush;
      RawVector dos_part(ncopy);
      RawVector hap_part(ncopy);
      std::memcpy(&dos_part[0], &gds_dos[k][0], (size_t)ncopy * sizeof(uint8_t));
      std::memcpy(&hap_part[0], &gds_hap[k][0], (size_t)ncopy * sizeof(uint8_t));
      dos_part.attr("dim") = IntegerVector::create(nsamples, since_flush);
      hap_part.attr("dim") = IntegerVector::create(nsamples, since_flush);
      dos_list[k] = dos_part;
      hap_list[k] = hap_part;
    }
    CharacterVector rchrom = r_chrom[Range(0, since_flush - 1)];
    CharacterVector rid    = r_id[Range(0, since_flush - 1)];
    CharacterVector rref   = r_ref[Range(0, since_flush - 1)];
    CharacterVector ralt   = r_alt[Range(0, since_flush - 1)];
    IntegerVector rpos     = r_pos[Range(0, since_flush - 1)];
    gds_append_fn(rchrom, rpos, rid, rref, ralt, dos_list, hap_list);
  }
}

// [[Rcpp::export]]
SEXP gds_to_txt_open_writer_cpp(std::string out_path, bool gz) {
  Writer* w = gz ? static_cast<Writer*>(new GzWriter(out_path))
                 : static_cast<Writer*>(new PlainWriter(out_path));
  Rcpp::XPtr<Writer> ptr(w, true);
  return ptr;
}

// [[Rcpp::export]]
void gds_to_txt_write_header_cpp(SEXP writer_ptr, CharacterVector sample_ids) {
  Rcpp::XPtr<Writer> writer(writer_ptr);
  std::string header = "CHROM\tPOS\tID\tREF\tALT";
  for (int i = 0; i < sample_ids.size(); ++i) {
    header += "\t" + std::string(CHAR(sample_ids[i]));
  }
  header += "\n";
  writer->write(header);
  writer->flush();
}

// [[Rcpp::export]]
void gds_to_txt_append_chunk_cpp(SEXP writer_ptr,
                                      CharacterVector chroms,
                                      IntegerVector pos,
                                      CharacterVector ids,
                                      CharacterVector refs,
                                      CharacterVector alts,
                                      NumericMatrix mat) {
  Rcpp::XPtr<Writer> writer(writer_ptr);
  const int nvars = mat.ncol();
  const int nsamples = mat.nrow();
  const double* m = mat.begin();
  std::string line;
  line.reserve(64 + nsamples * 2);
  for (int i = 0; i < nvars; ++i) {
    line.clear();
    line += std::string(CHAR(chroms[i])) + "\t" + std::to_string(pos[i]) + "\t" +
            std::string(CHAR(ids[i])) + "\t" + std::string(CHAR(refs[i])) + "\t" + std::string(CHAR(alts[i]));
    const double* col = m + (size_t)i * (size_t)nsamples;
    for (int j = 0; j < nsamples; ++j) {
      line.push_back('\t');
      line.push_back(static_cast<char>('0' + static_cast<int>(col[j])));
    }
    line.push_back('\n');
    writer->write(line);
  }
}

// [[Rcpp::export]]
void gds_to_txt_close_writer_cpp(SEXP writer_ptr) {
  Rcpp::XPtr<Writer> writer(writer_ptr);
  if (writer.get() != nullptr) {
    writer->flush();
    delete writer.get();
    R_ClearExternalPtr(writer_ptr);
  }
}

// [[Rcpp::export]]
void merge_task_results_cpp(std::string output_path,
                            std::string task_file_prefix,
                            std::string task_file_suffix,
                            int num_tasks) {
  auto writer = open_writer(output_path, ends_with(output_path, ".gz"));
  bool wrote_header = false;
  std::string line;

  int digits = 1;
  int tmp = num_tasks;
  while (tmp >= 10) {
    tmp /= 10;
    ++digits;
  }
  std::string task_file;
  task_file.reserve(task_file_prefix.size() + digits + task_file_suffix.size());

  for (int task_id = 1; task_id <= num_tasks; ++task_id) {
    task_file.assign(task_file_prefix);
    task_file.append(std::to_string(task_id));
    task_file.append(task_file_suffix);

    std::unique_ptr<LineReader> reader;
    try {
      reader = open_reader(task_file);
    } catch (...) {
      Rcpp::Rcout << "Warning: Task file not found: " << task_file << std::endl;
      continue;
    }

    bool first_line = true;
    while (reader->getline(line)) {
      if (first_line) {
        if (!wrote_header) {
          writer->write(line);
          writer->write("\n");
          wrote_header = true;
        }
        first_line = false;
        continue;
      }
      writer->write(line);
      writer->write("\n");
    }
  }

  writer->flush();
}

// ============================
// Linear regression helpers
// ============================

// Cholesky-first linear solver with QR fallback
// **Design:** Try Cholesky on X'X first (faster). If fails or ill-conditioned, fall back to pivoted QR.
// **Parameters:**
//   - use_qr_fallback: if false, R_out and perm_out are not returned when QR is used (beta still computed)
// **Returns:** 
//   - beta: solution vector (in original coordinate order)
//   - method: 'C' = Cholesky succeeded, 'Q' = pivoted QR used, 'F' = failed
//   - R_out: Upper Cholesky factor (Cholesky) or R factor (QR if use_qr_fallback=true) for SE computation
//   - perm_out: Column permutation indices (QR only, returned if use_qr_fallback=true)
static inline char solve_cholesky_qr_linear_eigen(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                   const Eigen::Ref<const Eigen::VectorXd>& y,
                                                   Eigen::VectorXd& beta,
                                                   Eigen::MatrixXd& XtX_work,
                                                   Eigen::VectorXd& Xty_work,
                                                   double tol_scale_rank = 1.0,
                                                   Eigen::MatrixXd* R_out = nullptr,
                                                   Eigen::VectorXi* perm_out = nullptr,
                                                   bool use_qr_fallback = true) {
  const int p = static_cast<int>(X.cols());
  beta = Eigen::VectorXd::Constant(p, NA_REAL);

  // Try Cholesky decomposition first using caller-provided work buffers
  XtX_work.setZero();
  XtX_work.template selfadjointView<Eigen::Lower>().rankUpdate(X.transpose());
  Xty_work.noalias() = X.transpose() * y;

  Eigen::LLT<Eigen::MatrixXd> llt(XtX_work);
  if (llt.info() == Eigen::Success) {
    // Validate Cholesky: check absolute floor and condition number
    Eigen::VectorXd L_diag = Eigen::MatrixXd(llt.matrixL()).diagonal().cwiseAbs();
    double min_diag = L_diag.minCoeff();
    double max_diag = L_diag.maxCoeff();
    
    if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
      // Cholesky succeeded and is well-conditioned
      beta = llt.solve(Xty_work);
      if (beta.allFinite()) {
        // Return R factor (upper Cholesky factor) for caller to compute SE via row-wise sums
        if (R_out) {
          *R_out = Eigen::MatrixXd(llt.matrixU());
        }
        // Note: no need to set perm_out for Cholesky (identity permutation is implicit)
        // Caller can assume identity permutation when method == 'C'
        return 'C';  // Cholesky succeeded
      }
    }
  }

  // Cholesky failed or ill-conditioned; try QR
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
  set_colpiv_qr_threshold(qr, X, tol_scale_rank);
  
  if (static_cast<int>(qr.rank()) < p) {
    return 'F';  // Failed: rank deficient
  }
  
  beta = qr.solve(y);
  if (!beta.allFinite()) {
    return 'F';  // Failed
  }
  
  // Store QR factors for SE calculation only if use_qr_fallback=true
  if (use_qr_fallback && R_out && perm_out) {
    *R_out = qr.matrixR().topLeftCorner(p, p).template triangularView<Eigen::Upper>();
    *perm_out = qr.colsPermutation().indices();
  }
  
  return 'Q';  // QR succeeded
}

// [[Rcpp::export]]
List tlstractor_linear_precompute(NumericVector y, NumericMatrix A,
                                  double tol_scale_rank = 1.0) {
  if (y.size() != A.nrow()) stop("y and A must have the same number of rows.");

  Eigen::Map<const Eigen::VectorXd> yv(y.begin(), y.size());
  Eigen::Map<const Eigen::MatrixXd> Am(A.begin(), A.nrow(), A.ncol());

  // Assume rankA > 0 (caller should ensure A has at least one non-zero column)
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Am);
  set_colpiv_qr_threshold(qr, Am, tol_scale_rank);

  const int n = static_cast<int>(Am.rows());
  int rankA = static_cast<int>(qr.rank());
  Eigen::MatrixXd A_reordered(n, rankA);
  Eigen::MatrixXd Qr(n, rankA);
  Eigen::MatrixXd R1(rankA, rankA);
  Eigen::VectorXd beta(rankA);
  Eigen::VectorXd eta = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd y_res = yv;
  Eigen::VectorXd Qty(rankA);
  Eigen::MatrixXd AtA_inv = Eigen::MatrixXd::Constant(rankA, rankA, NA_REAL);
  Eigen::MatrixXd inv_GammaAA = Eigen::MatrixXd::Constant(rankA, rankA, NA_REAL);
  Eigen::MatrixXd vcov_A_null = Eigen::MatrixXd::Constant(rankA, rankA, NA_REAL);

  if (rankA <= 0) {
    return List::create(
      Named("A") = wrap(A_reordered),
      Named("Qr") = wrap(Qr),
      Named("R1") = wrap(R1),
      Named("rank") = rankA,
      Named("beta") = wrap(beta),
      Named("eta") = wrap(eta),
      Named("y") = y,
      Named("y_res") = wrap(y_res),
      Named("Qty") = wrap(Qty),
      Named("AtA_inv") = wrap(AtA_inv),
      Named("inv_GammaAA") = wrap(inv_GammaAA),
      Named("vcov_A_null") = wrap(vcov_A_null),
      Named("R1_not_ok") = true
    );
  }

  Qr.setIdentity();
  qr.householderQ().setLength(rankA).applyThisOnTheLeft(Qr);
  R1 = qr.matrixR().topLeftCorner(rankA, rankA).template triangularView<Eigen::Upper>();

  Eigen::VectorXi pivot = qr.colsPermutation().indices();

  // Compute Qty = Qr^T * y, fitted null model eta = Qr * Qty, and residualized y: y_res = y - eta
  Qty = Qr.transpose() * yv;
  eta.noalias() = Qr * Qty;  // eta = A_reordered * beta in null model (projection onto col(A))
  y_res.noalias() -= eta;

  // Reorder A to match pivot order and truncate to rankA
  for (int j = 0; j < rankA; ++j) {
    A_reordered.col(j) = Am.col(pivot(j));
  }

  // Compute beta = R1^{-1} * Qty (coefficients from fitting y ~ A_reordered)
  beta = R1.template triangularView<Eigen::Upper>().solve(Qty);

  // Check R1 conditioning before inverse-based outputs
  Eigen::VectorXd R1_diag_abs = R1.diagonal().cwiseAbs();
  double R1_diag_min = R1_diag_abs.minCoeff();
  double R1_diag_max = R1_diag_abs.maxCoeff();
  double R1_diag_ratio = (R1_diag_min > 0.0) ? (R1_diag_max / R1_diag_min) : std::numeric_limits<double>::infinity();
  bool R1_not_ok = !(R1_diag_min > CHOLESKY_DIAG_FLOOR && R1_diag_ratio < CHOLESKY_COND_LIMIT);

  // Precompute (A'A)^{-1} = (R1'R1)^{-1} for efficient vcov computation (guarded by R1 check)
  if (!R1_not_ok) {
    AtA_inv.setIdentity(rankA, rankA);
    R1.transpose().triangularView<Eigen::Lower>().solveInPlace(AtA_inv);  // Solve R1^T * Y = I
    R1.triangularView<Eigen::Upper>().solveInPlace(AtA_inv);              // Solve R1 * X = Y
    inv_GammaAA.noalias() = static_cast<double>(n) * AtA_inv;  // inv((1/n) * A'A)
  }

  // Null-model vcov for reordered A: sigma2_null * (A'A)^{-1}
  const int df_resid_null = n - rankA;
  if (!R1_not_ok && df_resid_null > 0) {
    double sigma2_null = y_res.squaredNorm() / static_cast<double>(df_resid_null);
    if (std::isfinite(sigma2_null)) {
      vcov_A_null.noalias() = sigma2_null * AtA_inv;
    }
  }

  return List::create(
    Named("A") = wrap(A_reordered),    // n x rankA (reordered, rank-truncated)
    Named("Qr") = wrap(Qr),            // n x rankA
    Named("R1") = wrap(R1),            // rankA x rankA
    Named("rank") = rankA,
    Named("beta") = wrap(beta),        // rankA x 1, coefficients from y ~ A_reordered
    Named("eta") = wrap(eta),          // n x 1, A_reordered * beta in null model
    Named("y") = y,
    Named("y_res") = wrap(y_res),
    Named("Qty") = wrap(Qty),          // rankA x 1, used for efficient alpha computation
    Named("AtA_inv") = wrap(AtA_inv),  // rankA x rankA, (A'A)^{-1} for vcov computation
    Named("inv_GammaAA") = wrap(inv_GammaAA), // rankA x rankA, inv((1/n) * A'A)
    Named("vcov_A_null") = wrap(vcov_A_null),  // rankA x rankA null-model vcov for reordered A
    Named("R1_not_ok") = R1_not_ok
  );
}

struct LinearCoefAllWorkspaceEigen {
  Eigen::MatrixXd QtX;
  Eigen::MatrixXd X_res;
  Eigen::MatrixXd XtX_work;
  Eigen::VectorXd Xty_work;
  Eigen::VectorXd beta_x;
  Eigen::VectorXd Qty_adj;
  Eigen::VectorXd alpha;
  Eigen::VectorXd beta_all;
  bool used_qr_fallback = false;

  void resize(int n, int rankA, int pX) {
    QtX.resize(rankA, pX);
    X_res.resize(n, pX);
    XtX_work.resize(pX, pX);
    Xty_work.resize(pX);
    beta_x.resize(pX);
    Qty_adj.resize(rankA);
    alpha.resize(rankA);
    beta_all.resize(rankA + pX);
  }
};

struct LinearCoefXOnlyWorkspaceEigen {
  Eigen::MatrixXd XtX_work;
  Eigen::VectorXd Xty_work;
  Eigen::VectorXd beta;
  bool used_qr_fallback = false;

  void resize(int pX) {
    XtX_work.resize(pX, pX);
    Xty_work.resize(pX);
    beta.resize(pX);
  }
};

static inline bool tlstractor_linear_coef_all(
    const Eigen::Ref<const Eigen::MatrixXd>& Xm,
    const Eigen::Ref<const Eigen::MatrixXd>& Qr,
    const Eigen::Ref<const Eigen::MatrixXd>& R1,
    int rankA,
    const Eigen::Ref<const Eigen::VectorXd>& y_res,
    const Eigen::Ref<const Eigen::VectorXd>& Qty,
    double tol_scale_rank,
    LinearCoefAllWorkspaceEigen& ws) {
  const int pX = static_cast<int>(Xm.cols());

  ws.QtX.noalias() = Qr.transpose() * Xm;
  ws.X_res.noalias() = Xm;
  ws.X_res.noalias() -= Qr * ws.QtX;

  char method = solve_cholesky_qr_linear_eigen(
    ws.X_res,
    y_res,
    ws.beta_x,
    ws.XtX_work,
    ws.Xty_work,
    tol_scale_rank,
    nullptr,
    nullptr,
    true
  );
  if (method == 'F') {
    ws.beta_all.setConstant(NA_REAL);
    ws.used_qr_fallback = false;
    return false;
  }
  ws.used_qr_fallback = (method == 'Q');

  ws.Qty_adj.noalias() = Qty;
  ws.Qty_adj.noalias() -= ws.QtX * ws.beta_x;
  ws.alpha.noalias() = R1.template triangularView<Eigen::Upper>().solve(ws.Qty_adj);

  ws.beta_all.head(pX) = ws.beta_x;
  ws.beta_all.tail(rankA) = ws.alpha;
  return true;
}

static inline bool tlstractor_linear_coef_x_only(
    const Eigen::Ref<const Eigen::VectorXd>& yv,
    const Eigen::Ref<const Eigen::MatrixXd>& Xm,
    double tol_scale_rank,
    LinearCoefXOnlyWorkspaceEigen& ws) {
  char method = solve_cholesky_qr_linear_eigen(
    Xm,
    yv,
    ws.beta,
    ws.XtX_work,
    ws.Xty_work,
    tol_scale_rank,
    nullptr,
    nullptr,
    true
  );
  if (method == 'F') {
    ws.beta.setConstant(NA_REAL);
    ws.used_qr_fallback = false;
    return false;
  }
  ws.used_qr_fallback = (method == 'Q');
  return true;
}

struct LinearXStatsWorkspaceEigen {
  Eigen::MatrixXd QtX;
  Eigen::MatrixXd X_res;
  Eigen::MatrixXd XtX_work;
  Eigen::VectorXd Xty_work;
  Eigen::VectorXd y_work;
  Eigen::VectorXd beta;
  Eigen::MatrixXd R_x;
  Eigen::VectorXi perm_x;
  Eigen::VectorXd resid;
  Eigen::MatrixXd inv_R;
  Eigen::VectorXd se;
  Eigen::VectorXd se_decomp;
  Eigen::VectorXi perm_x_inv;
  Eigen::MatrixXd A_wald;
  Eigen::MatrixXd M_SS;
  Eigen::VectorXd beta_perm;

  void resize(int n, int rankA, int p, int df_wald) {
    QtX.resize(rankA, p);
    X_res.resize(n, p);
    XtX_work.resize(p, p);
    Xty_work.resize(p);
    y_work.resize(n);
    beta.resize(p);
    R_x.resize(p, p);
    perm_x.resize(p);
    resid.resize(n);
    inv_R.resize(p, p);
    se.resize(p);
    se_decomp.resize(p);
    perm_x_inv.resize(p);
    A_wald.resize(df_wald, p);
    M_SS.resize(df_wald, df_wald);
    beta_perm.resize(p);
  }
};

struct LinearXStatsResult {
  int df_resid;
  double wald;
  bool used_qr_fallback;
};

static inline LinearXStatsResult tlstractor_linear_x_stats_core(
    const Eigen::Ref<const Eigen::MatrixXd>& X_res,
    const Eigen::Ref<const Eigen::VectorXd>& y_work,
    int rankA,
    int df,
    double tol_scale_rank,
    bool use_qr_fallback,
    LinearXStatsWorkspaceEigen& ws) {
  const int n = static_cast<int>(X_res.rows());
  const int p = static_cast<int>(X_res.cols());

  ws.perm_x.setLinSpaced(p, 0, p - 1);
  char method = solve_cholesky_qr_linear_eigen(
    X_res,
    y_work,
    ws.beta,
    ws.XtX_work,
    ws.Xty_work,
    tol_scale_rank,
    &ws.R_x,
    &ws.perm_x,
    use_qr_fallback
  );

  if (method == 'F') {
    ws.beta.setConstant(NA_REAL);
    ws.se.setConstant(NA_REAL);
    return {NA_INTEGER, NA_REAL, false};
  }

  const int df_resid = n - rankA - p;
  if (df_resid <= 0) {
    ws.se.setConstant(NA_REAL);
    return {df_resid, NA_REAL, (method == 'Q')};
  }

  ws.resid.noalias() = y_work;
  ws.resid.noalias() -= X_res * ws.beta;
  const double sigma2 = ws.resid.squaredNorm() / static_cast<double>(df_resid);

  if (method == 'Q' && !use_qr_fallback) {
    ws.se.setConstant(NA_REAL);
    return {df_resid, NA_REAL, true};
  }

  ws.inv_R.setIdentity();
  ws.R_x.template triangularView<Eigen::Upper>().solveInPlace(ws.inv_R);
  ws.se.setZero();
  ws.se_decomp = (sigma2 * ws.inv_R.topRows(p).rowwise().squaredNorm().array()).sqrt().matrix();
  if (method == 'Q') {
    for (int j = 0; j < p; ++j) {
      ws.se(ws.perm_x(j)) = ws.se_decomp(j);
    }
  } else {
    ws.se = ws.se_decomp;
  }
  double wald_stat = NA_REAL;
  if (df == p) {
    if (method == 'Q') {
      ws.beta_perm.setZero();
      for (int j = 0; j < p; ++j) {
        ws.beta_perm(j) = ws.beta(ws.perm_x(j));
      }
      Eigen::VectorXd R_beta = ws.R_x.triangularView<Eigen::Upper>() * ws.beta_perm;
      wald_stat = R_beta.squaredNorm() / sigma2;
    } else {
      Eigen::VectorXd R_beta = ws.R_x.triangularView<Eigen::Upper>() * ws.beta;
      wald_stat = R_beta.squaredNorm() / sigma2;
    }
  } else if (df > 0 && df < p) {
    if (method == 'Q') {
      ws.perm_x_inv.setZero();
      for (int i = 0; i < p; ++i) {
        ws.perm_x_inv(ws.perm_x(i)) = i;
      }
      ws.A_wald.setZero();
      for (int j = 0; j < df; ++j) {
        int decomp_row = ws.perm_x_inv(j);
        ws.A_wald.row(j) = ws.inv_R.row(decomp_row);
      }
    } else {
      ws.A_wald.setZero();
      for (int j = 0; j < df; ++j) {
        ws.A_wald.row(j) = ws.inv_R.row(j);
      }
    }

    ws.M_SS = (ws.A_wald.topLeftCorner(df, p) * ws.A_wald.topLeftCorner(df, p).transpose()).eval();
    Eigen::VectorXd beta_S = ws.beta.head(df);

    Eigen::LLT<Eigen::MatrixXd> llt_M(ws.M_SS);
    if (llt_M.info() == Eigen::Success) {
      Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_M.matrixL()).diagonal().cwiseAbs();
      double min_diag = L_diag.minCoeff();
      double max_diag = L_diag.maxCoeff();
      if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
        Eigen::VectorXd x = llt_M.solve(beta_S);
        if (x.allFinite()) {
          double wald_raw = beta_S.dot(x) / sigma2;
          wald_stat = (wald_raw >= 0.0 && std::isfinite(wald_raw)) ? wald_raw : NA_REAL;
        }
      }
    }
  }

  return {df_resid, wald_stat, (method == 'Q')};
}

static inline LinearXStatsResult tlstractor_linear_x_stats(
    const Eigen::Ref<const Eigen::MatrixXd>& Xm,
    int df,
    const Eigen::Ref<const Eigen::MatrixXd>& Qr,
    const Eigen::Ref<const Eigen::VectorXd>& y_res,
    int rankA,
    double tol_scale_rank,
    bool use_qr_fallback,
    LinearXStatsWorkspaceEigen& ws) {
  ws.QtX.noalias() = Qr.transpose() * Xm;
  ws.X_res.noalias() = Xm;
  ws.X_res.noalias() -= Qr * ws.QtX;
  ws.y_work.noalias() = y_res;
  return tlstractor_linear_x_stats_core(ws.X_res, ws.y_work, rankA, df, tol_scale_rank, use_qr_fallback, ws);
}

static inline LinearXStatsResult tlstractor_linear_x_stats_only(
    const Eigen::Ref<const Eigen::MatrixXd>& Xm,
    int df,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    double tol_scale_rank,
    bool use_qr_fallback,
    LinearXStatsWorkspaceEigen& ws) {
  ws.X_res.noalias() = Xm;
  ws.y_work.noalias() = y;
  return tlstractor_linear_x_stats_core(ws.X_res, ws.y_work, 0, df, tol_scale_rank, use_qr_fallback, ws);
}

struct LinearCoefVcovAWorkspaceEigen {
  Eigen::VectorXd Qtx;
  Eigen::VectorXd x_res;
  Eigen::VectorXd Qty_adj;
  Eigen::VectorXd v1;
  Eigen::VectorXd resid;
  Eigen::VectorXd alpha;
  Eigen::MatrixXd vcov_alpha;
  double sigma2 = NA_REAL;
  double beta_x = NA_REAL;

  void resize(int n, int rankA) {
    Qtx.resize(rankA);
    x_res.resize(n);
    Qty_adj.resize(rankA);
    v1.resize(rankA);
    resid.resize(n);
    alpha.resize(rankA);
    vcov_alpha.resize(rankA, rankA);
  }
};

static inline bool tlstractor_linear_coef_vcov_A(
    const Eigen::Ref<const Eigen::VectorXd>& xv,
    const Eigen::Ref<const Eigen::MatrixXd>& Qr,
    const Eigen::Ref<const Eigen::MatrixXd>& R1,
    const Eigen::Ref<const Eigen::VectorXd>& y_res,
    const Eigen::Ref<const Eigen::VectorXd>& Qty,
    const Eigen::Ref<const Eigen::MatrixXd>& AtA_inv,
    int rankA,
    double beta_x,
    bool use_offset,
    LinearCoefVcovAWorkspaceEigen& ws) {
  const int n = static_cast<int>(y_res.size());
  ws.sigma2 = NA_REAL;
  ws.beta_x = beta_x;

  ws.Qtx.noalias() = Qr.transpose() * xv;
  ws.x_res.noalias() = xv;
  ws.x_res.noalias() -= Qr * ws.Qtx;

  if (!use_offset) {
    const double s_star = ws.x_res.squaredNorm();
    if (s_star < TOL_COLLINEARITY) {
      ws.alpha.setConstant(NA_REAL);
      ws.vcov_alpha.setConstant(NA_REAL);
      return false;
    }

    beta_x = ws.x_res.dot(y_res) / s_star;
    ws.beta_x = beta_x;
    ws.Qty_adj.noalias() = Qty;
    ws.Qty_adj.noalias() -= ws.Qtx * beta_x;
    ws.alpha = R1.triangularView<Eigen::Upper>().solve(ws.Qty_adj);

    ws.v1 = R1.triangularView<Eigen::Upper>().solve(ws.Qtx);
    ws.vcov_alpha.noalias() = AtA_inv;
    ws.vcov_alpha.noalias() += (ws.v1 * ws.v1.transpose()) / s_star;
  } else {
    ws.Qty_adj.noalias() = Qty;
    ws.Qty_adj.noalias() -= ws.Qtx * beta_x;
    ws.alpha = R1.triangularView<Eigen::Upper>().solve(ws.Qty_adj);
    ws.vcov_alpha.noalias() = AtA_inv;
  }

  ws.resid.noalias() = y_res;
  ws.resid.noalias() -= ws.x_res * beta_x;

  const int df_resid = n - rankA - (use_offset ? 0 : 1);
  double sigma2 = NA_REAL;
  if (df_resid > 0) {
    sigma2 = ws.resid.squaredNorm() / static_cast<double>(df_resid);
    ws.sigma2 = sigma2;
    ws.vcov_alpha *= sigma2;
  } else {
    ws.alpha.setConstant(NA_REAL);
    ws.vcov_alpha.setConstant(NA_REAL);
    return false;
  }

  return true;
}

// ============================
// tlstractor helpers
// ============================

static inline NumericMatrix subset_cols_cast_center_to_numeric(const SEXP mat_sexp,
                                                               const IntegerVector& idx,
                                                               const char* label) {
  const int nrow = Rf_nrows(mat_sexp);
  const int n_keep = idx.size();
  NumericMatrix out(nrow, n_keep);

  if (TYPEOF(mat_sexp) == INTSXP) {
    const int* in_ptr = INTEGER(mat_sexp);
    for (int j = 0; j < n_keep; ++j) {
      const int col_idx = idx[j] - 1;
      const int* in_col = in_ptr + static_cast<R_xlen_t>(col_idx) * nrow;
      double* out_col = out.begin() + static_cast<R_xlen_t>(j) * nrow;
      Eigen::Map<const Eigen::ArrayXi> in_col_map(in_col, nrow);
      Eigen::Map<Eigen::ArrayXd> out_col_map(out_col, nrow);
      out_col_map = in_col_map.cast<double>();
      out_col_map -= out_col_map.mean();
    }
    return out;
  }

  if (TYPEOF(mat_sexp) == REALSXP) {
    const double* in_ptr = REAL(mat_sexp);
    for (int j = 0; j < n_keep; ++j) {
      const int col_idx = idx[j] - 1;
      const double* in_col = in_ptr + static_cast<R_xlen_t>(col_idx) * nrow;
      double* out_col = out.begin() + static_cast<R_xlen_t>(j) * nrow;
      Eigen::Map<const Eigen::ArrayXd> in_col_map(in_col, nrow);
      Eigen::Map<Eigen::ArrayXd> out_col_map(out_col, nrow);
      out_col_map = in_col_map;
      out_col_map -= out_col_map.mean();
    }
    return out;
  }

  stop("Matrices in %s must be integer or numeric", label);
  return out;
}

static inline bool solve_spd_or_qr_eigen(const Eigen::Ref<const Eigen::MatrixXd>& M,
                                         const Eigen::Ref<const Eigen::VectorXd>& rhs,
                                         Eigen::VectorXd& sol,
                                         bool* used_qr_out) {
  *used_qr_out = false;
  Eigen::LLT<Eigen::MatrixXd> llt(M);
  if (llt.info() == Eigen::Success) {
    Eigen::VectorXd L_diag = Eigen::MatrixXd(llt.matrixL()).diagonal().cwiseAbs();
    const double min_diag = L_diag.minCoeff();
    const double max_diag = L_diag.maxCoeff();
    if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
      sol = llt.solve(rhs);
      if (sol.allFinite()) {
        return true;
      }
    }
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(M);
  set_colpiv_qr_threshold(qr, M, 1.0);
  if (static_cast<int>(qr.rank()) < static_cast<int>(M.cols())) {
    return false;
  }
  sol = qr.solve(rhs);
  *used_qr_out = true;
  return sol.allFinite();
}

static inline bool invert_spd_or_qr_eigen(const Eigen::Ref<const Eigen::MatrixXd>& M,
                                          Eigen::MatrixXd& Minv,
                                          bool* used_qr_out) {
  *used_qr_out = false;
  const int p = static_cast<int>(M.rows());

  Eigen::LLT<Eigen::MatrixXd> llt(M);
  if (llt.info() == Eigen::Success) {
    Eigen::VectorXd L_diag = Eigen::MatrixXd(llt.matrixL()).diagonal().cwiseAbs();
    const double min_diag = L_diag.minCoeff();
    const double max_diag = L_diag.maxCoeff();
    if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
      Minv.setIdentity(p, p);
      Eigen::MatrixXd L = Eigen::MatrixXd(llt.matrixL());
      L.template triangularView<Eigen::Lower>().solveInPlace(Minv);
      L.transpose().template triangularView<Eigen::Upper>().solveInPlace(Minv);
      Minv = Minv.template selfadjointView<Eigen::Lower>();
      return Minv.allFinite();
    }
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(M);
  set_colpiv_qr_threshold(qr, M, 1.0);
  if (static_cast<int>(qr.rank()) < p) {
    return false;
  }
  Minv.setIdentity(p, p);
  Minv = qr.solve(Minv);
  *used_qr_out = true;
  if (!Minv.allFinite()) {
    return false;
  }
  Minv = 0.5 * (Minv + Minv.transpose());
  return true;
}

struct GMMWorkspace {
  Eigen::MatrixXd X_full;
  Eigen::VectorXd Z;
  Eigen::VectorXd X_beta;
  Eigen::VectorXd XR_theta;
  Eigen::MatrixXd DX;
  Eigen::VectorXd DZ;
  Eigen::VectorXd DDA;
  Eigen::VectorXd delta12;
  Eigen::VectorXd gammaAZ;
  Eigen::VectorXd gammas;
  Eigen::MatrixXd inv_C;
  Eigen::MatrixXd Cn;
  Eigen::MatrixXd pseudo_X0;
  Eigen::MatrixXd ps_XtX;
  Eigen::VectorXd ps_y;
  Eigen::VectorXd ps_Xty;
  Eigen::VectorXd beta_full;
  Eigen::MatrixXd vcov_beta;
  Eigen::VectorXd solve_tmp;
  Eigen::VectorXd mu_prime;
  Eigen::VectorXd coeff_A;
  Eigen::MatrixXd XR_aug;
  Eigen::MatrixXd XR_weighted;
  Eigen::MatrixXd Gram_AA;
  Eigen::MatrixXd Gram_XR;
  Eigen::MatrixXd inv_Gram_AA;
  Eigen::MatrixXd inv_Gram_XR;

  void resize(int n, int p_full, int rankA, bool needs_softfail_buffers = false) {
    X_full.resize(n, p_full);
    Z.resize(n);
    X_beta.resize(n);
    XR_theta.resize(n);
    DX.resize(n, p_full);
    DZ.resize(n);
    DDA.resize(n);
    delta12.resize(p_full);
    gammaAZ.resize(rankA);
    gammas.resize(rankA);
    inv_C.resize(p_full + 1, p_full + 1);
    Cn.resize(p_full + 1, p_full + 1);
    pseudo_X0.resize(p_full + 1, p_full);
    ps_XtX.resize(p_full, p_full);
    ps_y.resize(p_full + 1);
    ps_Xty.resize(p_full);
    beta_full.resize(p_full);
    vcov_beta.resize(p_full, p_full);
    solve_tmp.resize(p_full, rankA);
    if (needs_softfail_buffers) {
      mu_prime.resize(n);
      coeff_A.resize(rankA);
      XR_aug.resize(n, rankA + 1);
      XR_weighted.resize(n, rankA + 1);
      Gram_AA.resize(rankA, rankA);
      Gram_XR.resize(rankA + 1, rankA + 1);
      inv_Gram_AA.resize(rankA, rankA);
      inv_Gram_XR.resize(rankA + 1, rankA + 1);
    }
  }
};

static inline bool delta_opt_gwas_linear_softpass(const Eigen::Ref<const Eigen::VectorXd>& y,
                                         const Eigen::Ref<const Eigen::VectorXd>& Z,
                                         const Eigen::Ref<const Eigen::MatrixXd>& X,
                                         const Eigen::Ref<const Eigen::MatrixXd>& A,
                                         const Eigen::Ref<const Eigen::VectorXd>& X_beta,
                                         const Eigen::Ref<const Eigen::VectorXd>& XR_theta,
                                         const Eigen::Ref<const Eigen::VectorXd>& A_thetaA,
                                         const Eigen::Ref<const Eigen::MatrixXd>& V_thetaA,
                                         double V_thetaZ,
                                         const Eigen::Ref<const Eigen::MatrixXd>& inv_GammaAA,
                                         GMMWorkspace& ws) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  const int rankA = static_cast<int>(A.cols());
  const double inv_n = 1.0 / static_cast<double>(n);

  ws.DX.noalias() = X;
  ws.DX.array().colwise() *= (X_beta - y).array();

  ws.DZ.noalias() = (Z.array() * (X_beta - XR_theta).array()).matrix();

  ws.inv_C.topLeftCorner(p, p).noalias() = inv_n * (ws.DX.transpose() * ws.DX);
  ws.delta12.noalias() = inv_n * (ws.DX.transpose() * ws.DZ);

  const double V_U2 = inv_n * ws.DZ.squaredNorm();
  const double GammaZZ = inv_n * Z.dot(Z);
  double Delta22 = V_U2 + (GammaZZ * GammaZZ) * (static_cast<double>(n) * V_thetaZ);

  if (rankA > 0) {
    ws.gammaAZ.noalias() = inv_n * (A.transpose() * Z);
    ws.gammas.noalias() = inv_GammaAA * ws.gammaAZ;
    ws.DDA.noalias() = A * ws.gammas;
    ws.DDA.array() *= (A_thetaA - y).array();

    ws.delta12.noalias() += inv_n * (ws.DX.transpose() * ws.DDA);

    const double Cov_U2theta = inv_n * ws.DZ.dot(ws.DDA);
    const double V_theta_term = static_cast<double>(n) *
      (ws.gammaAZ.transpose() * V_thetaA * ws.gammaAZ)(0, 0);
    Delta22 += V_theta_term + 2.0 * Cov_U2theta;
  }

  ws.inv_C.topRightCorner(p, 1) = ws.delta12;
  ws.inv_C.bottomLeftCorner(1, p) = ws.delta12.transpose();
  ws.inv_C(p, p) = Delta22;
  ws.inv_C = 0.5 * (ws.inv_C + ws.inv_C.transpose());

  return ws.inv_C.allFinite();
}

static inline bool delta_opt_gwas_logistic_softpass(const Eigen::Ref<const Eigen::VectorXd>& y,
                                           const Eigen::Ref<const Eigen::VectorXd>& Z,
                                           const Eigen::Ref<const Eigen::MatrixXd>& X,
                                           const Eigen::Ref<const Eigen::MatrixXd>& A,
                                           const Eigen::Ref<const Eigen::VectorXd>& mu_beta,
                                           const Eigen::Ref<const Eigen::VectorXd>& mu_theta,
                                           const Eigen::Ref<const Eigen::VectorXd>& mu_A_thetaA,
                                           const Eigen::Ref<const Eigen::MatrixXd>& V_thetaA,
                                           double V_thetaZ,
                                           const Eigen::Ref<const Eigen::MatrixXd>& inv_GammaAA,
                                           GMMWorkspace& ws) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  const double inv_n = 1.0 / static_cast<double>(n);

  Eigen::VectorXd& DZ = ws.DDA;
  Eigen::VectorXd& DZ2 = ws.DZ;
  Eigen::VectorXd& DDA = ws.XR_theta;

  ws.DX.noalias() = X;
  ws.DX.array().colwise() *= (mu_beta - y).array();

  DZ.noalias() = (Z.array() * (mu_beta - mu_theta).array()).matrix();
  DZ2.noalias() = (Z.array() * (mu_theta.array() * (1.0 - mu_theta.array())).array()).matrix();

  ws.inv_C.topLeftCorner(p, p).noalias() = inv_n * (ws.DX.transpose() * ws.DX);
  ws.delta12.noalias() = inv_n * (ws.DX.transpose() * DZ);

  const double V_U2 = inv_n * DZ.squaredNorm();
  const double GammaZZ = inv_n * DZ2.dot(Z);
  double Delta22 = V_U2 + (GammaZZ * GammaZZ) * (static_cast<double>(n) * V_thetaZ);

  ws.gammaAZ.noalias() = inv_n * (A.transpose() * DZ2);
  ws.gammas.noalias() = inv_GammaAA * ws.gammaAZ;
  DDA.noalias() = A * ws.gammas;
  DDA.array() *= (mu_A_thetaA - y).array();

  ws.delta12.noalias() += inv_n * (ws.DX.transpose() * DDA);

  const double Cov_U2theta = inv_n * DZ.dot(DDA);
  const double V_theta_term = static_cast<double>(n) *
    (ws.gammaAZ.transpose() * V_thetaA * ws.gammaAZ)(0, 0);
  Delta22 += V_theta_term + 2.0 * Cov_U2theta;

  ws.inv_C.topRightCorner(p, 1) = ws.delta12;
  ws.inv_C.bottomLeftCorner(1, p) = ws.delta12.transpose();
  ws.inv_C(p, p) = Delta22;
  ws.inv_C = 0.5 * (ws.inv_C + ws.inv_C.transpose());
  return ws.inv_C.allFinite();
}

static inline bool delta_opt_gwas_linear_softfail(const Eigen::Ref<const Eigen::VectorXd>& y,
                                                  const Eigen::Ref<const Eigen::VectorXd>& Z,
                                                  const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                  const Eigen::Ref<const Eigen::MatrixXd>& A,
                                                  const Eigen::Ref<const Eigen::VectorXd>& X_beta,
                                                  const Eigen::Ref<const Eigen::VectorXd>& XR_theta,
                                                  const Eigen::Ref<const Eigen::MatrixXd>& V_thetaA,
                                                  double V_thetaZ,
                                                  bool use_offset,
                                                  GMMWorkspace& ws,
                                                  bool* used_qr_out) {
  *used_qr_out = false;
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  const int rankA = static_cast<int>(A.cols());
  const double inv_n = 1.0 / static_cast<double>(n);

  ws.DX.noalias() = X;
  ws.DX.array().colwise() *= (X_beta - y).array();

  ws.DZ.noalias() = (Z.array() * (X_beta - XR_theta).array()).matrix();

  ws.inv_C.topLeftCorner(p, p).noalias() = inv_n * (ws.DX.transpose() * ws.DX);
  ws.delta12.noalias() = inv_n * (ws.DX.transpose() * ws.DZ);

  const double V_U2 = inv_n * ws.DZ.squaredNorm();
  const double GammaZZ = inv_n * Z.dot(Z);
  double Delta22 = V_U2 + (GammaZZ * GammaZZ) * (static_cast<double>(n) * V_thetaZ);

  if (rankA > 0) {
    ws.gammaAZ.noalias() = inv_n * (A.transpose() * Z);

    if (use_offset) {
      ws.Gram_AA.noalias() = inv_n * (A.transpose() * A);
      if (!invert_spd_or_qr_eigen(ws.Gram_AA, ws.inv_Gram_AA, used_qr_out)) {
        return false;
      }
      ws.gammas.noalias() = ws.inv_Gram_AA * ws.gammaAZ;
      ws.DDA.noalias() = A * ws.gammas;
      ws.DDA.array() *= (XR_theta - y).array();
    } else {
      ws.XR_aug.leftCols(rankA).noalias() = A;
      ws.XR_aug.col(rankA).noalias() = Z;

      ws.Gram_XR.noalias() = inv_n * (ws.XR_aug.transpose() * ws.XR_aug);
      if (!invert_spd_or_qr_eigen(ws.Gram_XR, ws.inv_Gram_XR, used_qr_out)) {
        return false;
      }
      ws.coeff_A.noalias() = ws.inv_Gram_XR.leftCols(rankA) * ws.gammaAZ;
      ws.DDA.noalias() = ws.XR_aug * ws.coeff_A;
      ws.DDA.array() *= (XR_theta - y).array();
    }

    ws.delta12.noalias() += inv_n * (ws.DX.transpose() * ws.DDA);

    const double Cov_U2theta = inv_n * ws.DZ.dot(ws.DDA);
    const double V_theta_term = static_cast<double>(n) *
      (ws.gammaAZ.transpose() * V_thetaA * ws.gammaAZ)(0, 0);
    Delta22 += V_theta_term + 2.0 * Cov_U2theta;
  }

  ws.inv_C.topRightCorner(p, 1) = ws.delta12;
  ws.inv_C.bottomLeftCorner(1, p) = ws.delta12.transpose();
  ws.inv_C(p, p) = Delta22;
  ws.inv_C = 0.5 * (ws.inv_C + ws.inv_C.transpose());
  return ws.inv_C.allFinite();
}

static inline bool delta_opt_gwas_logistic_softfail(const Eigen::Ref<const Eigen::VectorXd>& y,
                                                    const Eigen::Ref<const Eigen::VectorXd>& Z,
                                                    const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                    const Eigen::Ref<const Eigen::MatrixXd>& A,
                                                    const Eigen::Ref<const Eigen::VectorXd>& mu_beta,
                                                    const Eigen::Ref<const Eigen::VectorXd>& mu_XR,
                                                    const Eigen::Ref<const Eigen::MatrixXd>& V_thetaA,
                                                    double V_thetaZ,
                                                    bool use_offset,
                                                    GMMWorkspace& ws,
                                                    bool* used_qr_out) {
  *used_qr_out = false;
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  const int rankA = static_cast<int>(A.cols());
  const double inv_n = 1.0 / static_cast<double>(n);

  Eigen::VectorXd& DZ = ws.DDA;
  Eigen::VectorXd& DZ2 = ws.DZ;
  Eigen::VectorXd& DDA = ws.XR_theta;

  ws.mu_prime.noalias() = (mu_XR.array() * (1.0 - mu_XR.array())).matrix();

  ws.DX.noalias() = X;
  ws.DX.array().colwise() *= (mu_beta - y).array();

  DZ.noalias() = (Z.array() * (mu_beta - mu_XR).array()).matrix();
  DZ2.noalias() = (Z.array() * ws.mu_prime.array()).matrix();

  ws.inv_C.topLeftCorner(p, p).noalias() = inv_n * (ws.DX.transpose() * ws.DX);
  ws.delta12.noalias() = inv_n * (ws.DX.transpose() * DZ);

  const double V_U2 = inv_n * DZ.squaredNorm();
  const double GammaZZ = inv_n * DZ2.dot(Z);
  double Delta22 = V_U2 + (GammaZZ * GammaZZ) * (static_cast<double>(n) * V_thetaZ);

  ws.gammaAZ.noalias() = inv_n * (A.transpose() * DZ2);

  if (use_offset) {
    ws.XR_weighted.leftCols(rankA).noalias() = A;
    ws.XR_weighted.leftCols(rankA).array().colwise() *= ws.mu_prime.array();
    ws.Gram_AA.noalias() = inv_n * (A.transpose() * ws.XR_weighted.leftCols(rankA));
    if (!invert_spd_or_qr_eigen(ws.Gram_AA, ws.inv_Gram_AA, used_qr_out)) {
      return false;
    }
    ws.gammas.noalias() = ws.inv_Gram_AA * ws.gammaAZ;
    DDA.noalias() = A * ws.gammas;
    DDA.array() *= (mu_XR - y).array();
  } else {
    ws.XR_aug.leftCols(rankA).noalias() = A;
    ws.XR_aug.col(rankA).noalias() = Z;

    ws.XR_weighted.noalias() = ws.XR_aug;
    ws.XR_weighted.array().colwise() *= ws.mu_prime.array();

    ws.Gram_XR.noalias() = inv_n * (ws.XR_aug.transpose() * ws.XR_weighted);
    if (!invert_spd_or_qr_eigen(ws.Gram_XR, ws.inv_Gram_XR, used_qr_out)) {
      return false;
    }
    ws.coeff_A.noalias() = ws.inv_Gram_XR.leftCols(rankA) * ws.gammaAZ;
    DDA.noalias() = ws.XR_aug * ws.coeff_A;
    DDA.array() *= (mu_XR - y).array();
  }

  ws.delta12.noalias() += inv_n * (ws.DX.transpose() * DDA);

  const double Cov_U2theta = inv_n * DZ.dot(DDA);
  const double V_theta_term = static_cast<double>(n) *
    (ws.gammaAZ.transpose() * V_thetaA * ws.gammaAZ)(0, 0);
  Delta22 += V_theta_term + 2.0 * Cov_U2theta;

  ws.inv_C.topRightCorner(p, 1) = ws.delta12;
  ws.inv_C.bottomLeftCorner(1, p) = ws.delta12.transpose();
  ws.inv_C(p, p) = Delta22;
  ws.inv_C = 0.5 * (ws.inv_C + ws.inv_C.transpose());
  return ws.inv_C.allFinite();
}

static inline void direct_fast_build_pseudo_linear(const Eigen::Ref<const Eigen::MatrixXd>& Cn,
                                                   const Eigen::Ref<const Eigen::VectorXd>& y,
                                                   const Eigen::Ref<const Eigen::VectorXd>& Z,
                                                   const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                   const Eigen::Ref<const Eigen::VectorXd>& XR_theta,
                                                   GMMWorkspace& ws) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  const double scale_factor_beta = static_cast<double>(n) /
    std::sqrt(static_cast<double>(Cn.rows()));

  ws.pseudo_X0.topRows(p).noalias() = (X.transpose() * X) / scale_factor_beta;
  ws.pseudo_X0.row(p).noalias() = (Z.transpose() * X) / scale_factor_beta;

  ws.ps_XtX.noalias() = ws.pseudo_X0.transpose() * Cn * ws.pseudo_X0;

  ws.ps_y.head(p).noalias() = X.transpose() * y;
  ws.ps_y(p) = Z.dot(XR_theta);
  ws.ps_y /= scale_factor_beta;

  ws.ps_Xty.noalias() = ws.pseudo_X0.transpose() * Cn * ws.ps_y;
}

static inline void direct_fast_build_psXtX_linear(const Eigen::Ref<const Eigen::MatrixXd>& Cn,
                                                  const Eigen::Ref<const Eigen::VectorXd>& Z,
                                                  const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                  GMMWorkspace& ws) {
  const int n = static_cast<int>(X.rows());
  const int p = static_cast<int>(X.cols());
  const double scale_factor_beta = static_cast<double>(n) /
    std::sqrt(static_cast<double>(Cn.rows()));

  ws.pseudo_X0.topRows(p).noalias() = (X.transpose() * X) / scale_factor_beta;
  ws.pseudo_X0.row(p).noalias() = (Z.transpose() * X) / scale_factor_beta;
  ws.ps_XtX.noalias() = ws.pseudo_X0.transpose() * Cn * ws.pseudo_X0;
}

static inline void direct_fast_build_pseudo_logistic(const Eigen::Ref<const Eigen::MatrixXd>& Cn,
                                                     const Eigen::Ref<const Eigen::VectorXd>& y,
                                                     const Eigen::Ref<const Eigen::VectorXd>& Z,
                                                     const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                     const Eigen::Ref<const Eigen::VectorXd>& eta_beta,
                                                     const Eigen::Ref<const Eigen::VectorXd>& mu_beta,
                                                     const Eigen::Ref<const Eigen::VectorXd>& mu_theta,
                                                     GMMWorkspace& ws) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  const double scale_factor_beta = static_cast<double>(n) /
    std::sqrt(static_cast<double>(Cn.rows()));

  ws.DDA.noalias() = (mu_beta.array() * (1.0 - mu_beta.array())).matrix();
  ws.DDA.array() = ws.DDA.array().max(machine_eps());

  ws.DX.noalias() = X;
  ws.DX.array().colwise() *= ws.DDA.array();

  ws.pseudo_X0.topRows(p).noalias() = (X.transpose() * ws.DX) / scale_factor_beta;
  ws.pseudo_X0.row(p).noalias() = (Z.transpose() * ws.DX) / scale_factor_beta;
  ws.ps_XtX.noalias() = ws.pseudo_X0.transpose() * Cn * ws.pseudo_X0;

  ws.DZ.noalias() = (ws.DDA.array() * eta_beta.array() - mu_beta.array()).matrix();
  ws.ps_y.head(p).noalias() = X.transpose() * ws.DZ;
  ws.ps_y.head(p).noalias() += X.transpose() * y;

  ws.ps_y(p) = Z.dot(ws.DZ) + Z.dot(mu_theta);
  ws.ps_y /= scale_factor_beta;

  ws.ps_Xty.noalias() = ws.pseudo_X0.transpose() * Cn * ws.ps_y;
}

static inline void direct_fast_build_psXtX_logistic(const Eigen::Ref<const Eigen::MatrixXd>& Cn,
                                                    const Eigen::Ref<const Eigen::VectorXd>& Z,
                                                    const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                    const Eigen::Ref<const Eigen::VectorXd>& mu_beta,
                                                    GMMWorkspace& ws) {
  const int n = static_cast<int>(X.rows());
  const int p = static_cast<int>(X.cols());
  const double scale_factor_var = static_cast<double>(n);

  ws.DDA.noalias() = (mu_beta.array() * (1.0 - mu_beta.array())).matrix();
  ws.DDA.array() = ws.DDA.array().max(machine_eps());

  ws.DX.noalias() = X;
  ws.DX.array().colwise() *= ws.DDA.array();

  ws.pseudo_X0.topRows(p).noalias() = (X.transpose() * ws.DX) / scale_factor_var;
  ws.pseudo_X0.row(p).noalias() = (Z.transpose() * ws.DX) / scale_factor_var;
  ws.ps_XtX.noalias() = ws.pseudo_X0.transpose() * Cn * ws.pseudo_X0;
}

static inline bool solve_beta_vcov_from_pseudo_linear(const Eigen::Ref<const Eigen::MatrixXd>& ps_XtX,
                                                      const Eigen::Ref<const Eigen::VectorXd>& ps_Xty,
                                                      int n_samples,
                                                      int cn_rows,
                                                      GMMWorkspace& ws,
                                                      bool* used_qr_out) {
  *used_qr_out = false;
  const int p = static_cast<int>(ps_XtX.cols());
  const double vcov_scale = static_cast<double>(cn_rows) / static_cast<double>(n_samples);

  Eigen::LLT<Eigen::MatrixXd> llt(ps_XtX);
  if (llt.info() == Eigen::Success) {
    Eigen::VectorXd L_diag = Eigen::MatrixXd(llt.matrixL()).diagonal().cwiseAbs();
    const double min_diag = L_diag.minCoeff();
    const double max_diag = L_diag.maxCoeff();
    if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
    ws.beta_full.noalias() = llt.solve(ps_Xty);
    ws.vcov_beta.setIdentity(p, p);
    ws.vcov_beta = llt.solve(ws.vcov_beta);
    ws.vcov_beta *= vcov_scale;
    return ws.beta_full.allFinite() && ws.vcov_beta.allFinite();
    }
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(ps_XtX);
  set_colpiv_qr_threshold(qr, ps_XtX, 1.0);
  if (static_cast<int>(qr.rank()) < p) {
    return false;
  }
  *used_qr_out = true;

  ws.beta_full.noalias() = qr.solve(ps_Xty);
  if (!ws.beta_full.allFinite()) {
    return false;
  }

  ws.vcov_beta.setIdentity(p, p);
  ws.vcov_beta = qr.solve(ws.vcov_beta);
  if (!ws.vcov_beta.allFinite()) {
    return false;
  }
  ws.vcov_beta = 0.5 * (ws.vcov_beta + ws.vcov_beta.transpose());
  ws.vcov_beta *= vcov_scale;
  return true;
}

// ============================
// Linear regression for tlstractor
// ============================

// [[Rcpp::export]]
List fill_chunk_with_sumstats_softpass_linear(List dos_list, List hap_list, IntegerVector idx,
                                              NumericVector sumstats_beta, NumericVector sumstats_se,
                                              List precomp, List control, bool has_covar,
                                              bool cond_local, int num_ancs) {
  const int n_snps = idx.size();
  const int n_samples = Rf_nrows(dos_list[0]);
  const int n_laeff = cond_local ? num_ancs - 1 : 0;
  const int n_dos_nonref = num_ancs - 1;
  const int p_x = 1 + n_dos_nonref + n_laeff; // X = [Z, dos(anc1..), hap(...)]

  NumericMatrix beta(n_snps, num_ancs);
  NumericMatrix se(n_snps, num_ancs);
  NumericMatrix zstat(n_snps, num_ancs);
  std::fill(beta.begin(), beta.end(), NA_REAL);
  std::fill(se.begin(), se.end(), NA_REAL);
  std::fill(zstat.begin(), zstat.end(), NA_REAL);
  NumericVector wald(n_snps, NA_REAL);
  LogicalVector used_qr_fallback(n_snps, false);

  NumericMatrix laeff, lase, laz;
  if (cond_local) {
    laeff = NumericMatrix(n_snps, n_laeff);
    lase = NumericMatrix(n_snps, n_laeff);
    laz = NumericMatrix(n_snps, n_laeff);
    std::fill(laeff.begin(), laeff.end(), NA_REAL);
    std::fill(lase.begin(), lase.end(), NA_REAL);
    std::fill(laz.begin(), laz.end(), NA_REAL);
  }

  std::vector<NumericMatrix> dos_mats;
  dos_mats.reserve(n_dos_nonref);
  NumericMatrix z_mat = subset_cols_cast_center_to_numeric(dos_list[0], idx, "dos_list");
  Eigen::Map<Eigen::MatrixXd> Zm(z_mat.begin(), n_samples, n_snps);
  for (int k = 1; k < num_ancs; ++k) {
    dos_mats.emplace_back(subset_cols_cast_center_to_numeric(dos_list[k], idx, "dos_list"));
    Eigen::Map<const Eigen::MatrixXd> dmap(dos_mats.back().begin(), n_samples, n_snps);
    Zm.noalias() += dmap;
  }

  std::vector<NumericMatrix> hap_mats;
  if (cond_local) {
    hap_mats.reserve(n_laeff);
    for (int k = 0; k < n_laeff; ++k) {
      hap_mats.emplace_back(subset_cols_cast_center_to_numeric(hap_list[k], idx, "hap_list"));
    }
  }

  const bool refine_C = control.containsElementNamed("refine_C")
    ? as<bool>(control["refine_C"]) : false;

  NumericVector y_r = as<NumericVector>(precomp["y"]);
  Eigen::Map<const Eigen::VectorXd> y_map(y_r.begin(), y_r.size());

  NumericMatrix A_r;
  NumericVector eta_r;
  NumericMatrix inv_GammaAA_r;
  NumericMatrix vcov_A_null_r;
  NumericMatrix Qr_r;
  NumericMatrix R1_r;
  NumericVector yres_r;
  NumericVector Qty_r;
  int rankA = 0;

  if (has_covar) {
    Qr_r = as<NumericMatrix>(precomp["Qr"]);
    R1_r = as<NumericMatrix>(precomp["R1"]);
    yres_r = as<NumericVector>(precomp["y_res"]);
    Qty_r = as<NumericVector>(precomp["Qty"]);
    A_r = as<NumericMatrix>(precomp["A"]);
    eta_r = as<NumericVector>(precomp["eta"]);
    inv_GammaAA_r = as<NumericMatrix>(precomp["inv_GammaAA"]);
    vcov_A_null_r = as<NumericMatrix>(precomp["vcov_A_null"]);
    rankA = as<int>(precomp["rank"]);
  }

  const int p_full = p_x + rankA;
  GMMWorkspace gmm_ws;
  gmm_ws.resize(n_samples, p_full, rankA, false);

  std::optional<Eigen::Map<const Eigen::MatrixXd>> A_map;
  std::optional<Eigen::Map<const Eigen::VectorXd>> eta_map;
  std::optional<Eigen::Map<const Eigen::MatrixXd>> inv_GammaAA_map;
  std::optional<Eigen::Map<const Eigen::MatrixXd>> vcov_A_null_map;
  std::optional<Eigen::Map<const Eigen::MatrixXd>> Qr_map;
  std::optional<Eigen::Map<const Eigen::MatrixXd>> R1_map;
  std::optional<Eigen::Map<const Eigen::VectorXd>> yres_map;
  std::optional<Eigen::Map<const Eigen::VectorXd>> Qty_map;

  const Eigen::MatrixXd A_empty(n_samples, 0);
  const Eigen::MatrixXd inv_GammaAA_empty(0, 0);
  const Eigen::MatrixXd vcov_A_null_empty(0, 0);

  if (has_covar) {
    A_map.emplace(A_r.begin(), A_r.nrow(), A_r.ncol());
    eta_map.emplace(eta_r.begin(), eta_r.size());
    inv_GammaAA_map.emplace(inv_GammaAA_r.begin(), inv_GammaAA_r.nrow(), inv_GammaAA_r.ncol());
    vcov_A_null_map.emplace(vcov_A_null_r.begin(), vcov_A_null_r.nrow(), vcov_A_null_r.ncol());
    Qr_map.emplace(Qr_r.begin(), Qr_r.nrow(), Qr_r.ncol());
    R1_map.emplace(R1_r.begin(), R1_r.nrow(), R1_r.ncol());
    yres_map.emplace(yres_r.begin(), yres_r.size());
    Qty_map.emplace(Qty_r.begin(), Qty_r.size());
    gmm_ws.X_full.rightCols(rankA) = *A_map;
  }

  LinearCoefAllWorkspaceEigen coef_all_ws;
  LinearCoefXOnlyWorkspaceEigen coef_xonly_ws;

  if (has_covar) {
    coef_all_ws.resize(n_samples, rankA, p_x);
  } else {
    coef_xonly_ws.resize(p_x);
  }

  Eigen::MatrixXd T_dos = Eigen::MatrixXd::Zero(num_ancs, num_ancs);
  T_dos(0, 0) = 1.0;
  for (int k = 1; k < num_ancs; ++k) {
    T_dos(k, 0) = 1.0;
    T_dos(k, k) = 1.0;
  }

  auto fill_W_matrix = [&](int j) {
    gmm_ws.X_full.col(0).noalias() = Zm.col(j);
    for (int k = 1; k < num_ancs; ++k) {
      const NumericMatrix& dmat = dos_mats[k - 1];
      const double* src_col = dmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
      std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(k).data());
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const NumericMatrix& hmat = hap_mats[k];
        const int col_out = 1 + n_dos_nonref + k;
        const double* src_col = hmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
        std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(col_out).data());
      }
    }
  };

  for (int j = 0; j < n_snps; ++j) {
    bool used_qr_this_snp = false;
    const double curr_beta = sumstats_beta[j];
    const double curr_se = sumstats_se[j];
    const double V_thetaZ = curr_se * curr_se;
    if (!std::isfinite(V_thetaZ)) {
      continue;
    }
    
    fill_W_matrix(j);
    gmm_ws.Z.noalias() = Zm.col(j);

    bool init_ok = false;
    if (has_covar) {
      init_ok = tlstractor_linear_coef_all(
        gmm_ws.X_full.leftCols(p_x),
        *Qr_map,
        *R1_map,
        rankA,
        *yres_map,
        *Qty_map,
        1.0,
        coef_all_ws
      );
      if (init_ok) {
        gmm_ws.beta_full.noalias() = coef_all_ws.beta_all;
        used_qr_this_snp = coef_all_ws.used_qr_fallback;
      }
    } else {
      init_ok = tlstractor_linear_coef_x_only(y_map, gmm_ws.X_full.leftCols(p_x), 1.0, coef_xonly_ws);
      if (init_ok) {
        gmm_ws.beta_full.noalias() = coef_xonly_ws.beta;
        used_qr_this_snp = coef_xonly_ws.used_qr_fallback;
      }
    }

    if (!init_ok || !gmm_ws.beta_full.allFinite()) {
      continue;
    }

    gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;
    if (has_covar) {
      gmm_ws.XR_theta.noalias() = *eta_map;
    } else {
      gmm_ws.XR_theta.setZero();
    }
    gmm_ws.XR_theta.noalias() += gmm_ws.Z * curr_beta;

    if (has_covar) {
      if (!delta_opt_gwas_linear_softpass(y_map, gmm_ws.Z, gmm_ws.X_full,
                                 *A_map,
                                 gmm_ws.X_beta, gmm_ws.XR_theta,
                                 *eta_map,
                                 *vcov_A_null_map,
                                 V_thetaZ,
                                 *inv_GammaAA_map,
                                 gmm_ws)) {
        continue;
      }
    } else {
      if (!delta_opt_gwas_linear_softpass(y_map, gmm_ws.Z, gmm_ws.X_full,
                                 A_empty,
                                 gmm_ws.X_beta, gmm_ws.XR_theta,
                                 y_map,
                                 vcov_A_null_empty,
                                 V_thetaZ,
                                 inv_GammaAA_empty,
                                 gmm_ws)) {
        continue;
      }
    }
    bool used_qr_tmp = false;
    if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

    direct_fast_build_pseudo_linear(gmm_ws.Cn, y_map, gmm_ws.Z, gmm_ws.X_full, gmm_ws.XR_theta, gmm_ws);
    if (!gmm_ws.ps_XtX.allFinite() || !gmm_ws.ps_Xty.allFinite()) {
      continue;
    }

    if (!refine_C) {
      used_qr_tmp = false;
      if (!solve_beta_vcov_from_pseudo_linear(gmm_ws.ps_XtX, gmm_ws.ps_Xty,
                                              n_samples, gmm_ws.Cn.rows(), gmm_ws, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
    } else {
      used_qr_tmp = false;
      if (!solve_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.ps_Xty, gmm_ws.beta_full, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

      gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;
      if (has_covar) {
        if (!delta_opt_gwas_linear_softpass(y_map, gmm_ws.Z, gmm_ws.X_full,
                                   *A_map,
                                   gmm_ws.X_beta, gmm_ws.XR_theta,
                                   *eta_map,
                                   *vcov_A_null_map,
                                   V_thetaZ,
                                   *inv_GammaAA_map,
                                   gmm_ws)) {
          continue;
        }
      } else {
        if (!delta_opt_gwas_linear_softpass(y_map, gmm_ws.Z, gmm_ws.X_full,
                                   A_empty,
                                   gmm_ws.X_beta, gmm_ws.XR_theta,
                                   y_map,
                                   vcov_A_null_empty,
                                   V_thetaZ,
                                   inv_GammaAA_empty,
                                   gmm_ws)) {
          continue;
        }
      }
      used_qr_tmp = false;
      if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

      direct_fast_build_psXtX_linear(gmm_ws.Cn, gmm_ws.Z, gmm_ws.X_full, gmm_ws);
      used_qr_tmp = false;
      if (!invert_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.vcov_beta, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
      gmm_ws.vcov_beta *= static_cast<double>(gmm_ws.Cn.rows()) / static_cast<double>(n_samples);
      if (!gmm_ws.vcov_beta.allFinite()) {
        continue;
      }
    }

    Eigen::VectorXd beta_dos_param = gmm_ws.beta_full.head(num_ancs);
    Eigen::MatrixXd cov_dos_param = gmm_ws.vcov_beta.topLeftCorner(num_ancs, num_ancs);
    Eigen::VectorXd beta_dos = T_dos * beta_dos_param;
    Eigen::MatrixXd cov_dos = T_dos * cov_dos_param * T_dos.transpose();
    cov_dos = 0.5 * (cov_dos + cov_dos.transpose());

    for (int k = 0; k < num_ancs; ++k) {
      beta(j, k) = beta_dos(k);
      const double vkk = cov_dos(k, k);
      if (std::isfinite(vkk) && vkk > 0.0) {
        se(j, k) = std::sqrt(vkk);
        const double z = beta(j, k) / se(j, k);
        zstat(j, k) = std::isfinite(z) ? z : NA_REAL;
      }
    }

    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const int col = 1 + n_dos_nonref + k;
        laeff(j, k) = gmm_ws.beta_full(col);
        const double vkk = gmm_ws.vcov_beta(col, col);
        if (std::isfinite(vkk) && vkk > 0.0) {
          lase(j, k) = std::sqrt(vkk);
          const double z = laeff(j, k) / lase(j, k);
          laz(j, k) = std::isfinite(z) ? z : NA_REAL;
        }
      }
    }

    Eigen::VectorXd wald_tmp(num_ancs);
    bool used_qr_wald = false;
    if (solve_spd_or_qr_eigen(cov_dos, beta_dos, wald_tmp, &used_qr_wald)) {
      used_qr_this_snp = used_qr_this_snp || used_qr_wald;
      const double wald_stat = beta_dos.dot(wald_tmp);
      if (std::isfinite(wald_stat) && wald_stat >= 0.0) {
        wald[j] = wald_stat;
      }
    }
    used_qr_fallback[j] = used_qr_this_snp;
  }

  List out = List::create(
    Named("beta") = beta,
    Named("se") = se,
    Named("z") = zstat,
    Named("wald") = wald,
    Named("used_qr_fallback") = used_qr_fallback
  );
  if (cond_local) {
    out["laeff"] = laeff;
    out["lase"] = lase;
    out["laz"] = laz;
  }
  return out;
}

// [[Rcpp::export]]
List fill_chunk_with_sumstats_softfail_linear(List dos_list, List hap_list, IntegerVector idx,
                                              NumericVector sumstats_beta, NumericVector sumstats_se,
                                              List precomp, List control, bool has_covar,
                                              bool cond_local, int num_ancs) {
  const int n_snps = idx.size();
  const int n_samples = Rf_nrows(dos_list[0]);
  const int n_laeff = cond_local ? num_ancs - 1 : 0;
  const int n_dos_nonref = num_ancs - 1;
  const int p_x = 1 + n_dos_nonref + n_laeff;

  NumericMatrix beta(n_snps, num_ancs);
  NumericMatrix se(n_snps, num_ancs);
  NumericMatrix zstat(n_snps, num_ancs);
  std::fill(beta.begin(), beta.end(), NA_REAL);
  std::fill(se.begin(), se.end(), NA_REAL);
  std::fill(zstat.begin(), zstat.end(), NA_REAL);
  NumericVector wald(n_snps, NA_REAL);
  LogicalVector used_qr_fallback(n_snps, false);

  NumericMatrix laeff, lase, laz;
  if (cond_local) {
    laeff = NumericMatrix(n_snps, n_laeff);
    lase = NumericMatrix(n_snps, n_laeff);
    laz = NumericMatrix(n_snps, n_laeff);
    std::fill(laeff.begin(), laeff.end(), NA_REAL);
    std::fill(lase.begin(), lase.end(), NA_REAL);
    std::fill(laz.begin(), laz.end(), NA_REAL);
  }

  std::vector<NumericMatrix> dos_mats;
  dos_mats.reserve(n_dos_nonref);
  NumericMatrix z_mat = subset_cols_cast_center_to_numeric(dos_list[0], idx, "dos_list");
  Eigen::Map<Eigen::MatrixXd> Zm(z_mat.begin(), n_samples, n_snps);
  for (int k = 1; k < num_ancs; ++k) {
    dos_mats.emplace_back(subset_cols_cast_center_to_numeric(dos_list[k], idx, "dos_list"));
    Eigen::Map<const Eigen::MatrixXd> dmap(dos_mats.back().begin(), n_samples, n_snps);
    Zm.noalias() += dmap;
  }

  std::vector<NumericMatrix> hap_mats;
  if (cond_local) {
    hap_mats.reserve(n_laeff);
    for (int k = 0; k < n_laeff; ++k) {
      hap_mats.emplace_back(subset_cols_cast_center_to_numeric(hap_list[k], idx, "hap_list"));
    }
  }

  const bool refine_C = control.containsElementNamed("refine_C")
    ? as<bool>(control["refine_C"]) : false;
  const bool use_offset = control.containsElementNamed("use_offset")
    ? as<bool>(control["use_offset"]) : true;

  NumericVector y_r = as<NumericVector>(precomp["y"]);
  Eigen::Map<const Eigen::VectorXd> y_map(y_r.begin(), y_r.size());

  NumericMatrix A_r;
  NumericMatrix AtA_inv_r;
  NumericMatrix Qr_r;
  NumericMatrix R1_r;
  NumericVector yres_r;
  NumericVector Qty_r;
  int rankA = 0;

  if (has_covar) {
    Qr_r = as<NumericMatrix>(precomp["Qr"]);
    R1_r = as<NumericMatrix>(precomp["R1"]);
    yres_r = as<NumericVector>(precomp["y_res"]);
    Qty_r = as<NumericVector>(precomp["Qty"]);
    A_r = as<NumericMatrix>(precomp["A"]);
    AtA_inv_r = as<NumericMatrix>(precomp["AtA_inv"]);
    rankA = as<int>(precomp["rank"]);
  }

  const int p_full = p_x + rankA;
  GMMWorkspace gmm_ws;
  gmm_ws.resize(n_samples, p_full, rankA, true);

  std::optional<Eigen::Map<const Eigen::MatrixXd>> A_map;
  std::optional<Eigen::Map<const Eigen::MatrixXd>> AtA_inv_map;
  std::optional<Eigen::Map<const Eigen::MatrixXd>> Qr_map;
  std::optional<Eigen::Map<const Eigen::MatrixXd>> R1_map;
  std::optional<Eigen::Map<const Eigen::VectorXd>> yres_map;
  std::optional<Eigen::Map<const Eigen::VectorXd>> Qty_map;

  const Eigen::MatrixXd A_empty(n_samples, 0);
  const Eigen::MatrixXd inv_GammaAA_empty(0, 0);
  const Eigen::MatrixXd vcov_A_empty(0, 0);

  if (has_covar) {
    A_map.emplace(A_r.begin(), A_r.nrow(), A_r.ncol());
    AtA_inv_map.emplace(AtA_inv_r.begin(), AtA_inv_r.nrow(), AtA_inv_r.ncol());
    Qr_map.emplace(Qr_r.begin(), Qr_r.nrow(), Qr_r.ncol());
    R1_map.emplace(R1_r.begin(), R1_r.nrow(), R1_r.ncol());
    yres_map.emplace(yres_r.begin(), yres_r.size());
    Qty_map.emplace(Qty_r.begin(), Qty_r.size());
    gmm_ws.X_full.rightCols(rankA) = *A_map;
  }

  LinearCoefAllWorkspaceEigen coef_all_ws;
  LinearCoefXOnlyWorkspaceEigen coef_xonly_ws;
  LinearCoefVcovAWorkspaceEigen coef_A_ws;

  if (has_covar) {
    coef_all_ws.resize(n_samples, rankA, p_x);
    coef_A_ws.resize(n_samples, rankA);
  } else {
    coef_xonly_ws.resize(p_x);
  }

  Eigen::MatrixXd T_dos = Eigen::MatrixXd::Zero(num_ancs, num_ancs);
  T_dos(0, 0) = 1.0;
  for (int k = 1; k < num_ancs; ++k) {
    T_dos(k, 0) = 1.0;
    T_dos(k, k) = 1.0;
  }

  auto fill_W_matrix = [&](int j) {
    gmm_ws.X_full.col(0).noalias() = Zm.col(j);
    for (int k = 1; k < num_ancs; ++k) {
      const NumericMatrix& dmat = dos_mats[k - 1];
      const double* src_col = dmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
      std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(k).data());
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const NumericMatrix& hmat = hap_mats[k];
        const int col_out = 1 + n_dos_nonref + k;
        const double* src_col = hmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
        std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(col_out).data());
      }
    }
  };

  for (int j = 0; j < n_snps; ++j) {
    bool used_qr_this_snp = false;
    const double curr_beta = sumstats_beta[j];
    const double curr_se = sumstats_se[j];
    const double V_thetaZ = curr_se * curr_se;
    if (!std::isfinite(V_thetaZ)) {
      continue;
    }

    fill_W_matrix(j);
    gmm_ws.Z.noalias() = Zm.col(j);

    bool init_ok = false;
    if (has_covar) {
      init_ok = tlstractor_linear_coef_all(
        gmm_ws.X_full.leftCols(p_x),
        *Qr_map,
        *R1_map,
        rankA,
        *yres_map,
        *Qty_map,
        1.0,
        coef_all_ws
      );
      if (init_ok) {
        gmm_ws.beta_full.noalias() = coef_all_ws.beta_all;
        used_qr_this_snp = coef_all_ws.used_qr_fallback;
      }
    } else {
      init_ok = tlstractor_linear_coef_x_only(y_map, gmm_ws.X_full.leftCols(p_x), 1.0, coef_xonly_ws);
      if (init_ok) {
        gmm_ws.beta_full.noalias() = coef_xonly_ws.beta;
        used_qr_this_snp = coef_xonly_ws.used_qr_fallback;
      }
    }

    if (!init_ok || !gmm_ws.beta_full.allFinite()) {
      continue;
    }

    gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;

    if (has_covar) {
      if (!tlstractor_linear_coef_vcov_A(
            gmm_ws.Z,
            *Qr_map,
            *R1_map,
            *yres_map,
            *Qty_map,
            *AtA_inv_map,
            rankA,
            curr_beta,
            use_offset,
            coef_A_ws)) {
        continue;
      }

      gmm_ws.XR_theta.noalias() = (*A_map) * coef_A_ws.alpha;
      gmm_ws.XR_theta.noalias() += gmm_ws.Z * curr_beta;

      bool used_qr_delta = false;
      if (!delta_opt_gwas_linear_softfail(y_map, gmm_ws.Z, gmm_ws.X_full,
                                          *A_map,
                                          gmm_ws.X_beta,
                                          gmm_ws.XR_theta,
                                          coef_A_ws.vcov_alpha,
                                          V_thetaZ,
                                          use_offset,
                                          gmm_ws,
                                          &used_qr_delta)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_delta;
    } else {
      gmm_ws.XR_theta.noalias() = gmm_ws.Z * curr_beta;
      if (!delta_opt_gwas_linear_softpass(y_map, gmm_ws.Z, gmm_ws.X_full,
                                 A_empty,
                                 gmm_ws.X_beta, gmm_ws.XR_theta,
                                 y_map,
                                 vcov_A_empty,
                                 V_thetaZ,
                                 inv_GammaAA_empty,
                                 gmm_ws)) {
        continue;
      }
    }

    bool used_qr_tmp = false;
    if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

    direct_fast_build_pseudo_linear(gmm_ws.Cn, y_map, gmm_ws.Z, gmm_ws.X_full, gmm_ws.XR_theta, gmm_ws);
    if (!gmm_ws.ps_XtX.allFinite() || !gmm_ws.ps_Xty.allFinite()) {
      continue;
    }

    if (!refine_C) {
      used_qr_tmp = false;
      if (!solve_beta_vcov_from_pseudo_linear(gmm_ws.ps_XtX, gmm_ws.ps_Xty,
                                              n_samples, gmm_ws.Cn.rows(), gmm_ws, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
    } else {
      used_qr_tmp = false;
      if (!solve_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.ps_Xty, gmm_ws.beta_full, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

      gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;

      if (has_covar) {
        bool used_qr_delta = false;
        if (!delta_opt_gwas_linear_softfail(y_map, gmm_ws.Z, gmm_ws.X_full,
                                            *A_map,
                                            gmm_ws.X_beta,
                                            gmm_ws.XR_theta,
                                            coef_A_ws.vcov_alpha,
                                            V_thetaZ,
                                            use_offset,
                                            gmm_ws,
                                            &used_qr_delta)) {
          continue;
        }
        used_qr_this_snp = used_qr_this_snp || used_qr_delta;
      } else {
        if (!delta_opt_gwas_linear_softpass(y_map, gmm_ws.Z, gmm_ws.X_full,
                                   A_empty,
                                   gmm_ws.X_beta, gmm_ws.XR_theta,
                                   y_map,
                                   vcov_A_empty,
                                   V_thetaZ,
                                   inv_GammaAA_empty,
                                   gmm_ws)) {
          continue;
        }
      }

      used_qr_tmp = false;
      if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

      direct_fast_build_psXtX_linear(gmm_ws.Cn, gmm_ws.Z, gmm_ws.X_full, gmm_ws);
      used_qr_tmp = false;
      if (!invert_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.vcov_beta, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
      gmm_ws.vcov_beta *= static_cast<double>(gmm_ws.Cn.rows()) / static_cast<double>(n_samples);
      if (!gmm_ws.vcov_beta.allFinite()) {
        continue;
      }
    }

    Eigen::VectorXd beta_dos_param = gmm_ws.beta_full.head(num_ancs);
    Eigen::MatrixXd cov_dos_param = gmm_ws.vcov_beta.topLeftCorner(num_ancs, num_ancs);
    Eigen::VectorXd beta_dos = T_dos * beta_dos_param;
    Eigen::MatrixXd cov_dos = T_dos * cov_dos_param * T_dos.transpose();
    cov_dos = 0.5 * (cov_dos + cov_dos.transpose());

    for (int k = 0; k < num_ancs; ++k) {
      beta(j, k) = beta_dos(k);
      const double vkk = cov_dos(k, k);
      if (std::isfinite(vkk) && vkk > 0.0) {
        se(j, k) = std::sqrt(vkk);
        const double z = beta(j, k) / se(j, k);
        zstat(j, k) = std::isfinite(z) ? z : NA_REAL;
      }
    }

    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const int col = 1 + n_dos_nonref + k;
        laeff(j, k) = gmm_ws.beta_full(col);
        const double vkk = gmm_ws.vcov_beta(col, col);
        if (std::isfinite(vkk) && vkk > 0.0) {
          lase(j, k) = std::sqrt(vkk);
          const double z = laeff(j, k) / lase(j, k);
          laz(j, k) = std::isfinite(z) ? z : NA_REAL;
        }
      }
    }

    Eigen::VectorXd wald_tmp(num_ancs);
    bool used_qr_wald = false;
    if (solve_spd_or_qr_eigen(cov_dos, beta_dos, wald_tmp, &used_qr_wald)) {
      used_qr_this_snp = used_qr_this_snp || used_qr_wald;
      const double wald_stat = beta_dos.dot(wald_tmp);
      if (std::isfinite(wald_stat) && wald_stat >= 0.0) {
        wald[j] = wald_stat;
      }
    }
    used_qr_fallback[j] = used_qr_this_snp;
  }

  List out = List::create(
    Named("beta") = beta,
    Named("se") = se,
    Named("z") = zstat,
    Named("wald") = wald,
    Named("used_qr_fallback") = used_qr_fallback
  );
  if (cond_local) {
    out["laeff"] = laeff;
    out["lase"] = lase;
    out["laz"] = laz;
  }
  return out;
}

// [[Rcpp::export]]
List fill_chunk_without_sumstats_linear(List dos_list, List hap_list, IntegerVector idx,
                                        List precomp, List control, bool has_covar, bool cond_local, int num_ancs) {
  const int n_snps = idx.size();
  const int n_samples = Rf_nrows(dos_list[0]);
  const int n_laeff = cond_local ? num_ancs - 1 : 0;
  const int p_total = num_ancs + n_laeff;

  NumericMatrix beta(n_snps, num_ancs);
  NumericMatrix se(n_snps, num_ancs);
  NumericMatrix tstat(n_snps, num_ancs);
  std::fill(beta.begin(), beta.end(), NA_REAL);
  std::fill(se.begin(), se.end(), NA_REAL);
  std::fill(tstat.begin(), tstat.end(), NA_REAL);
  NumericVector wald(n_snps, NA_REAL);
  LogicalVector used_qr_fallback(n_snps, false);
  IntegerVector df_resid(n_snps, NA_INTEGER);

  NumericMatrix laeff, lase, lat;
  if (cond_local) {
    laeff = NumericMatrix(n_snps, n_laeff);
    lase = NumericMatrix(n_snps, n_laeff);
    lat = NumericMatrix(n_snps, n_laeff);
    std::fill(laeff.begin(), laeff.end(), NA_REAL);
    std::fill(lase.begin(), lase.end(), NA_REAL);
    std::fill(lat.begin(), lat.end(), NA_REAL);
  }

  std::vector<NumericMatrix> dos_mats;
  dos_mats.reserve(num_ancs);
  for (int k = 0; k < num_ancs; ++k) {
    dos_mats.emplace_back(subset_cols_cast_center_to_numeric(dos_list[k], idx, "dos_list"));
  }

  std::vector<NumericMatrix> hap_mats;
  if (cond_local) {
    hap_mats.reserve(num_ancs);
    for (int k = 0; k < num_ancs; ++k) {
      hap_mats.emplace_back(subset_cols_cast_center_to_numeric(hap_list[k], idx, "hap_list"));
    }
  }

  const bool use_qr_fallback = control.containsElementNamed("use_qr_fallback")
    ? as<bool>(control["use_qr_fallback"]) : true;
  const int rankA = has_covar ? as<int>(precomp["rank"]) : 0;

  NumericMatrix Qr_r;
  NumericVector yres_r;
  NumericVector y_r;
  if (has_covar) {
    Qr_r = as<NumericMatrix>(precomp["Qr"]);
    yres_r = as<NumericVector>(precomp["y_res"]);
  } else {
    y_r = as<NumericVector>(precomp["y"]);
  }

  NumericMatrix X(n_samples, p_total);
  Eigen::Map<Eigen::MatrixXd> Xm(X.begin(), n_samples, p_total);

  std::optional<Eigen::Map<const Eigen::MatrixXd>> Qr_map;
  std::optional<Eigen::Map<const Eigen::VectorXd>> yres_map;
  std::optional<Eigen::Map<const Eigen::VectorXd>> y_map;
  if (has_covar) {
    Qr_map.emplace(Qr_r.begin(), Qr_r.nrow(), Qr_r.ncol());
    yres_map.emplace(yres_r.begin(), yres_r.size());
  } else {
    y_map.emplace(y_r.begin(), y_r.size());
  }

  const int p = p_total;
  const int df_wald = num_ancs;
  LinearXStatsWorkspaceEigen linear_workspace;
  linear_workspace.resize(n_samples, rankA, p, df_wald);

  for (int j = 0; j < n_snps; ++j) {
    for (int k = 0; k < num_ancs; ++k) {
      const NumericMatrix& dmat = dos_mats[k];
      const double* src_col = dmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
      double* dst_col = X.begin() + static_cast<R_xlen_t>(k) * n_samples;
      std::copy(src_col, src_col + n_samples, dst_col);
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const NumericMatrix& hmat = hap_mats[k];
        const int col_out = num_ancs + k;
        const double* src_col = hmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
        double* dst_col = X.begin() + static_cast<R_xlen_t>(col_out) * n_samples;
        std::copy(src_col, src_col + n_samples, dst_col);
      }
    }

    LinearXStatsResult stats;
    if (has_covar) {
      stats = tlstractor_linear_x_stats(
        Xm,
        num_ancs,
        *Qr_map,
        *yres_map,
        rankA,
        1.0,
        use_qr_fallback,
        linear_workspace
      );
    } else {
      stats = tlstractor_linear_x_stats_only(
        Xm,
        num_ancs,
        *y_map,
        1.0,
        use_qr_fallback,
        linear_workspace
      );
    }

    df_resid[j] = stats.df_resid;
    used_qr_fallback[j] = stats.used_qr_fallback;
    for (int k = 0; k < num_ancs; ++k) {
      beta(j, k) = linear_workspace.beta(k);
      se(j, k) = linear_workspace.se(k);
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const int col_in = num_ancs + k;
        laeff(j, k) = linear_workspace.beta(col_in);
        lase(j, k) = linear_workspace.se(col_in);
      }
    }
    if (std::isfinite(stats.wald) && stats.wald >= 0.0) {
      wald[j] = stats.wald;
    }
  }

  Eigen::Map<const Eigen::ArrayXXd> beta_arr(beta.begin(), n_snps, num_ancs);
  Eigen::Map<const Eigen::ArrayXXd> se_arr(se.begin(), n_snps, num_ancs);
  Eigen::Map<Eigen::ArrayXXd> t_arr(tstat.begin(), n_snps, num_ancs);
  t_arr = beta_arr / se_arr;

  if (cond_local) {
    Eigen::Map<const Eigen::ArrayXXd> laeff_arr(laeff.begin(), n_snps, n_laeff);
    Eigen::Map<const Eigen::ArrayXXd> lase_arr(lase.begin(), n_snps, n_laeff);
    Eigen::Map<Eigen::ArrayXXd> lat_arr(lat.begin(), n_snps, n_laeff);
    lat_arr = laeff_arr / lase_arr;
  }

  List out = List::create(
    Named("beta") = beta,
    Named("se") = se,
    Named("t") = tstat,
    Named("wald") = wald,
    Named("used_qr_fallback") = used_qr_fallback,
    Named("df_resid") = df_resid
  );
  if (cond_local) {
    out["laeff"] = laeff;
    out["lase"] = lase;
    out["lat"] = lat;
  }
  
  return out;
}

// ============================
// Logistic regression helpers
// ============================

static inline void logistic_mu_stable_eigen(const Eigen::Ref<const Eigen::VectorXd>& eta,
                                            Eigen::Ref<Eigen::VectorXd> mu) {
  for (Eigen::Index i = 0; i < eta.size(); ++i) {
    const double x = eta[i];
    double v;
    if (x >= 0.0) {
      const double ex = std::exp(-x);
      v = 1.0 / (1.0 + ex);
    } else {
      const double ex = std::exp(x);
      v = ex / (1.0 + ex);
    }
    mu[i] = v;
  }
}

// Numerically stable binomial deviance using softplus trick
// deviance = 2 * sum(max(0,eta) + log(1+exp(-|eta|)) - y*eta)
// This is more efficient and equally stable as the previous log-sum-exp approach
static inline double binom_deviance_eigen(const Eigen::Ref<const Eigen::VectorXd>& y,
                                          const Eigen::Ref<const Eigen::VectorXd>& eta) {
  double dev = 0.0;
  for (Eigen::Index i = 0; i < y.size(); ++i) {
    const double yi = y[i];
    const double etai = eta[i];
    // Softplus formula: max(0,eta) + log(1+exp(-|eta|)) = log(1+exp(eta))
    const double term = std::max(0.0, etai) + std::log1p(std::exp(-std::abs(etai))) - yi * etai;
    dev += 2.0 * term;
  }
  return dev;
}

// WLS solver: tries Cholesky first, falls back to pivoted QR if Cholesky fails or is ill-conditioned
// Pre-allocate working vectors/matrices in caller to avoid repeated allocations in loops.
static inline void wls_cholesky_qr_step_eigen(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                              const Eigen::Ref<const Eigen::VectorXd>& z,
                                              const Eigen::Ref<const Eigen::VectorXd>& W,
                                              Eigen::Ref<Eigen::VectorXd> sqrtW_work,
                                              Eigen::Ref<Eigen::MatrixXd> Xw_work,
                                              Eigen::Ref<Eigen::VectorXd> zw_work,
                                              Eigen::Ref<Eigen::MatrixXd> XtWX_work,
                                              Eigen::Ref<Eigen::VectorXd> XtWz_work,
                                              Eigen::Ref<Eigen::VectorXd> beta_out,
                                              int* rank_out,
                                              bool* used_qr,
                                              double tol_scale_rank = 1.0,
                                              Eigen::MatrixXd* R_out = nullptr,
                                              Eigen::VectorXi* perm_out = nullptr) {
  const int p = static_cast<int>(X.cols());
  
  // Reuse pre-allocated work vectors
  sqrtW_work.noalias() = W.array().sqrt().matrix();
  Xw_work.noalias() = X;
  Xw_work.array().colwise() *= sqrtW_work.array();
  zw_work.noalias() = (z.array() * sqrtW_work.array()).matrix();
  
  // Compute X'WX using rankUpdate: XtWX = Xw' * Xw (more efficient via self-adjoint rank update)
  XtWX_work.topLeftCorner(p, p).setZero();
  XtWX_work.topLeftCorner(p, p)
      .template selfadjointView<Eigen::Lower>()
      .rankUpdate(Xw_work.leftCols(p).transpose());
  
  // Compute X'Wz directly
  XtWz_work.head(p).noalias() = Xw_work.leftCols(p).transpose() * zw_work;
  
  Eigen::LLT<Eigen::MatrixXd> llt(XtWX_work.topLeftCorner(p, p));
  if (llt.info() == Eigen::Success) {
    // Check condition number via diagonal ratio
    Eigen::VectorXd L_diag = Eigen::MatrixXd(llt.matrixL()).diagonal().cwiseAbs();
    double min_diag = L_diag.minCoeff();
    double max_diag = L_diag.maxCoeff();
    if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
      beta_out.noalias() = llt.solve(XtWz_work.head(p));
      if (beta_out.allFinite()) {
        *rank_out = p;
        *used_qr = false;
        if (R_out) {
          *R_out = Eigen::MatrixXd(llt.matrixU());
        }
        return;
      }
    }
  }
  
  // Fallback to pivoted QR using work vectors
  *used_qr = true;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Xw_work.leftCols(p));
  set_colpiv_qr_threshold(qr, Xw_work.leftCols(p), tol_scale_rank);
  *rank_out = static_cast<int>(qr.rank());
  beta_out.noalias() = qr.solve(zw_work);
  if (R_out && perm_out) {
    *R_out = qr.matrixR().topLeftCorner(p, p).template triangularView<Eigen::Upper>();
    *perm_out = qr.colsPermutation().indices();
  }
}

// Hybrid Schur WLS step for IRLS with separated blocks: [X, A].
// Uses direct Schur complement: eliminate A block first, solve for beta_x, then recover beta_a.
// If L_cached is provided (non-null), uses it as the Cholesky factor of M.
// Returns {success, beta_x, beta_a, rank_deficient}.
struct SchurStepResult {
  bool success;
  Eigen::VectorXd beta_x;
  Eigen::VectorXd beta_a;
  bool rank_deficient;
  bool used_qr_for_solve;
};

static inline SchurStepResult wls_schur_step_eigen(const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                   const Eigen::Ref<const Eigen::MatrixXd>& A,
                                                   const Eigen::Ref<const Eigen::VectorXd>& z,
                                                   const Eigen::Ref<const Eigen::VectorXd>& W,
                                                   Eigen::Ref<Eigen::VectorXd> sqrtW_work,
                                                   Eigen::Ref<Eigen::MatrixXd> Aw_work,
                                                   Eigen::Ref<Eigen::MatrixXd> Xw_work,
                                                   Eigen::Ref<Eigen::VectorXd> zw_work,
                                                   Eigen::Ref<Eigen::MatrixXd> M_work,
                                                   Eigen::Ref<Eigen::MatrixXd> B_work,
                                                   Eigen::Ref<Eigen::VectorXd> c_work,
                                                   Eigen::Ref<Eigen::MatrixXd> MinvB_work,
                                                   Eigen::Ref<Eigen::VectorXd> Minvc_work,
                                                   Eigen::Ref<Eigen::MatrixXd> S_work,
                                                   Eigen::Ref<Eigen::VectorXd> r_work,
                                                   const Eigen::Ref<const Eigen::MatrixXd>* L_cached = nullptr,
                                                   bool L_cached_not_ok = false,
                                                   double tol_scale_rank = 1.0,
                                                   Eigen::MatrixXd* R_s_out = nullptr,
                                                   Eigen::VectorXi* perm_s_out = nullptr) {
  SchurStepResult res;
  res.success = false;
  res.rank_deficient = false;
  res.used_qr_for_solve = false;

  const int pX = static_cast<int>(X.cols());
  const int rankA = static_cast<int>(A.cols());

  // Form weighted matrices: Aw = sqrt(W) * A, etc. (reuse pre-allocated work matrices)
  sqrtW_work.noalias() = W.array().sqrt().matrix();
  Aw_work.noalias() = A;
  Aw_work.array().colwise() *= sqrtW_work.array();
  Xw_work.noalias() = X;
  Xw_work.array().colwise() *= sqrtW_work.array();
  zw_work.noalias() = (z.array() * sqrtW_work.array()).matrix();
  
  // Step 1: Get M = A'WA
  // Try Cholesky on M with proper thresholds
  // If L_cached is provided and valid, reuse it directly for solves
  Eigen::LLT<Eigen::MatrixXd> llt_M;
  const Eigen::Ref<const Eigen::MatrixXd>* L_ptr = nullptr;
  bool use_cholesky_M = false;
  
  if (L_cached != nullptr) {
    L_ptr = L_cached;
    use_cholesky_M = true;
  } else if (!L_cached_not_ok) {
    // Use rankUpdate for efficiency: only computes lower triangle
    M_work.topLeftCorner(rankA, rankA).setZero();
    M_work.topLeftCorner(rankA, rankA)
      .template selfadjointView<Eigen::Lower>()
      .rankUpdate(Aw_work.topLeftCorner(Aw_work.rows(), rankA).transpose());

    llt_M.compute(M_work);
    if (llt_M.info() == Eigen::Success) {
      Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_M.matrixL()).diagonal().cwiseAbs();
      double min_diag = L_diag.minCoeff();
      double max_diag = L_diag.maxCoeff();
      if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
        use_cholesky_M = true;
      }
    }
  }
  
  if (!use_cholesky_M) {
    // Fallback: pivoted QR on Aw to check rank
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_A(Aw_work.topLeftCorner(Aw_work.rows(), rankA));
    set_colpiv_qr_threshold(qr_A, Aw_work.topLeftCorner(Aw_work.rows(), rankA), tol_scale_rank);
    int rA = static_cast<int>(qr_A.rank());
    if (rA < rankA) {
      res.rank_deficient = true;
      return res;
    }
    // Rank is full but Cholesky failed -> signal to use fallback in caller
    return res;
  }
  
  // Step 2: Form B = Aw' * Xw and c = Aw' * zw
  B_work.noalias() = Aw_work.topLeftCorner(Aw_work.rows(), rankA).transpose() * Xw_work.leftCols(pX);  // rankA × pX
  c_work.noalias() = Aw_work.topLeftCorner(Aw_work.rows(), rankA).transpose() * zw_work;   // rankA × 1
  
  // Step 3: Solve M \ B and M \ c using Cholesky (does NOT form full inv(M))
  auto solve_M_mat = [&](const Eigen::MatrixXd& rhs) -> Eigen::MatrixXd {
    if (L_ptr) {
      Eigen::MatrixXd y = L_ptr->triangularView<Eigen::Lower>().solve(rhs);
      Eigen::MatrixXd x = L_ptr->transpose().triangularView<Eigen::Upper>().solve(y);
      return x;
    }
    Eigen::MatrixXd x = llt_M.solve(rhs);
    return x;
  };
  auto solve_M_vec = [&](const Eigen::VectorXd& rhs) -> Eigen::VectorXd {
    if (L_ptr) {
      Eigen::VectorXd y = L_ptr->triangularView<Eigen::Lower>().solve(rhs);
      Eigen::VectorXd x = L_ptr->transpose().triangularView<Eigen::Upper>().solve(y);
      return x;
    }
    Eigen::VectorXd x = llt_M.solve(rhs);
    return x;
  };

  MinvB_work.noalias() = solve_M_mat(B_work);    // Solves M * MinvB = B
  Minvc_work.noalias() = solve_M_vec(c_work);    // Solves M * Minvc = c
  
  // Step 4: Form Schur complement system
  // S = X'WX - B'(M^(-1))B  (computed via matrix products, not explicit inverse)
  // Use rankUpdate for first term (only computes lower triangle), then subtract second term
  S_work.topLeftCorner(pX, pX).setZero();
  S_work.topLeftCorner(pX, pX)
    .template selfadjointView<Eigen::Lower>()
    .rankUpdate(Xw_work.leftCols(pX).transpose());
  S_work.noalias() -= B_work.transpose() * MinvB_work;  // Subtract B'M^{-1}B (symmetric, full computation)
  S_work = S_work.template selfadjointView<Eigen::Lower>();
  
  // r = X'Wz - B'(M^(-1))c
  r_work.noalias() = Xw_work.leftCols(pX).transpose() * zw_work - B_work.transpose() * Minvc_work;  // pX × 1
  
  // Step 5: Solve S * beta_x = r with scale-aware Cholesky thresholds
  Eigen::LLT<Eigen::MatrixXd> llt_S(S_work);
  Eigen::VectorXd beta_x;
  bool use_cholesky_S = false;
  
  if (llt_S.info() == Eigen::Success) {
    Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_S.matrixL()).diagonal().cwiseAbs();
    double min_diag = L_diag.minCoeff();
    double max_diag = L_diag.maxCoeff();
    if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
      beta_x = llt_S.solve(r_work);
      use_cholesky_S = true;
      res.used_qr_for_solve = false;
      if (R_s_out) {
        *R_s_out = Eigen::MatrixXd(llt_S.matrixU());
      }
    }
  }
  
  if (!use_cholesky_S) {
    // Cholesky failed or ill-conditioned -> try pivoted QR to check rank
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_S(S_work);
    set_colpiv_qr_threshold(qr_S, S_work, tol_scale_rank);
    
    int rS = static_cast<int>(qr_S.rank());
    if (rS < pX) {
      res.rank_deficient = true;
      return res;
    }
    res.used_qr_for_solve = true;
    if (R_s_out && perm_s_out) {
      *R_s_out = qr_S.matrixR().topLeftCorner(pX, pX).template triangularView<Eigen::Upper>();
      *perm_s_out = qr_S.colsPermutation().indices();
    }
    beta_x = qr_S.solve(r_work);
  }
  
  // Step 6: Recover beta_a = M^(-1) * (c - B * beta_x) via solving, not forming inverse
  Eigen::VectorXd beta_a = solve_M_vec(c_work - B_work * beta_x);
  
  res.success = true;
  res.beta_x = beta_x;
  res.beta_a = beta_a;
  return res;
}

// [[Rcpp::export]]
List tlstractor_logistic_precompute(NumericVector y, NumericMatrix A,
                                    int maxit = 25, double tol = 1e-8,
                                    int max_step_halving = 10,
                                    double tol_scale_rank = 1.0) {
  if (y.size() != A.nrow()) stop("y and A must have the same number of rows.");

  Eigen::Map<const Eigen::VectorXd> yv(y.begin(), y.size());
  Eigen::Map<const Eigen::MatrixXd> Am(A.begin(), A.nrow(), A.ncol());

  // Check rank of A via pivoted QR (same as linear)
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Am);
  set_colpiv_qr_threshold(qr, Am, tol_scale_rank);
  int rankA = static_cast<int>(qr.rank());
  Eigen::VectorXi pivot = qr.colsPermutation().indices();
  
  // Reorder A according to pivot, keep only rankA columns
  Eigen::MatrixXd A_reordered(Am.rows(), rankA);
  for (int j = 0; j < rankA; ++j) {
    A_reordered.col(j) = Am.col(pivot(j));
  }

  Eigen::VectorXd beta = Eigen::VectorXd::Zero(rankA);
  Eigen::VectorXd eta = Eigen::VectorXd::Zero(yv.size());

  Eigen::VectorXd mu(yv.size());
  logistic_mu_stable_eigen(eta, mu);
  double dev = binom_deviance_eigen(yv, eta);

  bool converged = false;
  bool has_error = false;
  std::string error_msg = "";
  int it = 0;
  
  // Preallocate working vectors to avoid repeated allocations in loop
  Eigen::VectorXd W(yv.size());
  Eigen::VectorXd z(yv.size());
  Eigen::VectorXd eta_new(yv.size());
  Eigen::VectorXd mu_new(yv.size());
  Eigen::VectorXd beta_new(rankA);
  
  // Preallocate WLS work vectors/matrices for wls_cholesky_qr_step_eigen
  Eigen::VectorXd sqrtW_work(yv.size());
  Eigen::MatrixXd Xw_work(yv.size(), rankA);
  Eigen::VectorXd zw_work(yv.size());
  Eigen::MatrixXd XtWX_work(rankA, rankA);
  Eigen::VectorXd XtWz_work(rankA);

  for (; it < maxit; ++it) {
    // Calculate W and Z
    W.noalias() = (mu.array() * (1.0 - mu.array())).matrix();
    z.noalias() = eta + (yv - mu).cwiseQuotient(W);

    // Solve WLS: tries Cholesky first, falls back to QR
    int rank_wls = 0;
    bool used_qr = false;
    wls_cholesky_qr_step_eigen(A_reordered, z, W, sqrtW_work, Xw_work, zw_work, 
                               XtWX_work, XtWz_work, beta_new, &rank_wls, &used_qr, tol_scale_rank);
    
    // Check for rank deficiency (signals near-separation)
    if (rank_wls < rankA) {
      has_error = true;
      error_msg = "rank_deficient_in_irls";
      break;
    }

    // Compute new eta, mu, deviance
    eta_new.noalias() = A_reordered * beta_new;
    logistic_mu_stable_eigen(eta_new, mu_new);
    double dev_new = binom_deviance_eigen(yv, eta_new);

    // Check mu exactly 0 or 1 (vectorized)
    bool has_exact_boundary = (mu_new.array() <= 0.0).any() || (mu_new.array() >= 1.0).any();
    
    // Check for invalid deviance or worse deviance or exact boundary -> try step-halving
    bool need_halving = (!std::isfinite(dev_new)) || (dev_new > dev) || has_exact_boundary;
    
    if (need_halving) {
      int halving = 0;
      while (halving < max_step_halving) {
        beta_new.noalias() = 0.5 * (beta_new + beta);
        eta_new.noalias() = A_reordered * beta_new;
        logistic_mu_stable_eigen(eta_new, mu_new);
        dev_new = binom_deviance_eigen(yv, eta_new);
        
        // Check if improved (vectorized)
        has_exact_boundary = (mu_new.array() <= 0.0).any() || (mu_new.array() >= 1.0).any();
        
        if (std::isfinite(dev_new) && dev_new <= dev && !has_exact_boundary) {
          break;  // Step-halving succeeded
        }
        ++halving;
      }
      
      // After max halvings, check if still invalid
      if (!std::isfinite(dev_new) || has_exact_boundary) {
        has_error = true;
        error_msg = has_exact_boundary ? "mu_at_boundary_after_halving" : "infinite_deviance_after_halving";
        break;
      }
    }

    // Update and check convergence
    double max_diff = (beta_new - beta).cwiseAbs().maxCoeff();
    double rel_dev = std::abs(dev - dev_new) / (0.1 + std::abs(dev));

    beta.noalias() = beta_new;
    eta.noalias() = eta_new;
    dev = dev_new;
    mu.noalias() = mu_new;

    if (rel_dev < tol && max_diff < tol) {
      converged = true;
      break;
    }
  }

  // Compute L_AtWA at convergence for reuse in first subsequent Schur iteration
  // Also cache W and z for warm-start in first per-SNP iteration
  Eigen::MatrixXd AtWA = Eigen::MatrixXd::Zero(rankA, rankA);
  Eigen::MatrixXd L_AtWA = Eigen::MatrixXd::Zero(rankA, rankA);
  Eigen::MatrixXd vcov_A_null = Eigen::MatrixXd::Constant(rankA, rankA, NA_REAL);
  Eigen::MatrixXd inv_GammaAA = Eigen::MatrixXd::Constant(rankA, rankA, NA_REAL);
  bool L_AtWA_not_ok = true;
  W.setZero();
  z.setZero();
  
  if (converged) {
    W.noalias() = (mu.array() * (1.0 - mu.array())).matrix();
    z.noalias() = eta + (yv - mu).cwiseQuotient(W);
    
    Eigen::VectorXd sqrtW = W.array().sqrt().matrix();
    Eigen::MatrixXd Aw = A_reordered.array().colwise() * sqrtW.array();
    // Use rankUpdate for efficiency: only computes lower triangle
    AtWA.template selfadjointView<Eigen::Lower>().rankUpdate(Aw.transpose());
    AtWA = AtWA.template selfadjointView<Eigen::Lower>();
    
    // Compute Cholesky factor for caching
    Eigen::LLT<Eigen::MatrixXd> llt_AtWA(AtWA);
    if (llt_AtWA.info() == Eigen::Success) {
      Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_AtWA.matrixL()).diagonal().cwiseAbs();
      double min_diag = L_diag.minCoeff();
      double max_diag = L_diag.maxCoeff();
      if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
        L_AtWA = llt_AtWA.matrixL();
        L_AtWA_not_ok = false;

        // Null-model vcov for reordered A in logistic model: (A'WA)^{-1}
        vcov_A_null.setIdentity(rankA, rankA);
        L_AtWA.template triangularView<Eigen::Lower>().solveInPlace(vcov_A_null);
        L_AtWA.transpose().template triangularView<Eigen::Upper>().solveInPlace(vcov_A_null);
        vcov_A_null = vcov_A_null.template selfadjointView<Eigen::Lower>();
        inv_GammaAA.noalias() = static_cast<double>(Am.rows()) * vcov_A_null; // inv((1/n) * A'WA)
      }
    }
  }

  return List::create(
    Named("beta") = wrap(beta),             // rankA coefficients (reordered)
    Named("eta") = wrap(eta),
    Named("mu") = wrap(mu),
    Named("y") = y,
    Named("A") = wrap(A_reordered),         // n x rankA (reordered, rank-truncated)
    Named("rank") = rankA,
    Named("converged") = converged,
    Named("iter") = converged ? (it + 1) : maxit,
    Named("deviance") = dev,
    Named("error") = has_error,
    Named("error_msg") = error_msg,
    Named("L_AtWA") = wrap(L_AtWA),         // Cholesky factor of AtWA for reuse
    Named("L_AtWA_not_ok") = L_AtWA_not_ok,         // Whether L_AtWA is valid for reuse
    Named("vcov_A_null") = wrap(vcov_A_null),       // rankA x rankA null-model vcov for reordered A
    Named("inv_GammaAA") = wrap(inv_GammaAA),       // rankA x rankA, inv((1/n) * A'WA)
    Named("W") = wrap(W),                   // Weight vector at convergence for warm-start
    Named("z") = wrap(z)                    // Pseudo-response at convergence for warm-start
  );
}

struct LogisticCoefAllWorkspaceEigen {
  Eigen::VectorXd beta;
  Eigen::VectorXd eta;
  Eigen::VectorXd mu;
  Eigen::VectorXd W;
  Eigen::VectorXd z;
  Eigen::VectorXd eta_new;
  Eigen::VectorXd mu_new;
  Eigen::VectorXd beta_new;
  Eigen::MatrixXd Z;

  Eigen::VectorXd sqrtW_schur;
  Eigen::MatrixXd Aw_schur;
  Eigen::MatrixXd Xw_schur;
  Eigen::VectorXd zw_schur;
  Eigen::MatrixXd M_schur;
  Eigen::MatrixXd B_schur;
  Eigen::VectorXd c_schur;
  Eigen::MatrixXd MinvB_schur;
  Eigen::VectorXd Minvc_schur;
  Eigen::MatrixXd S_schur;
  Eigen::VectorXd r_schur;

  Eigen::VectorXd sqrtW_wls;
  Eigen::MatrixXd Xw_wls;
  Eigen::VectorXd zw_wls;
  Eigen::MatrixXd XtWX_wls;
  Eigen::VectorXd XtWz_wls;

  void resize(int n, int rankA, int pX, int p_full) {
    beta.resize(p_full);
    eta.resize(n);
    mu.resize(n);
    W.resize(n);
    z.resize(n);
    eta_new.resize(n);
    mu_new.resize(n);
    beta_new.resize(p_full);
    Z.resize(n, p_full);

    sqrtW_schur.resize(n);
    Aw_schur.resize(n, rankA);
    Xw_schur.resize(n, pX);
    zw_schur.resize(n);
    M_schur.resize(rankA, rankA);
    B_schur.resize(rankA, pX);
    c_schur.resize(rankA);
    MinvB_schur.resize(rankA, pX);
    Minvc_schur.resize(rankA);
    S_schur.resize(pX, pX);
    r_schur.resize(pX);

    sqrtW_wls.resize(n);
    Xw_wls.resize(n, p_full);
    zw_wls.resize(n);
    XtWX_wls.resize(p_full, p_full);
    XtWz_wls.resize(p_full);
  }
};

struct LogisticCoefAllResult {
  bool converged;
  int iter;
  double deviance;
  bool error;
  std::string error_msg;
  bool used_qr_fallback;
};

static inline LogisticCoefAllResult tlstractor_logistic_coef_all(
    const Eigen::Ref<const Eigen::MatrixXd>& Xm,
    const Eigen::Ref<const Eigen::MatrixXd>& Am,
    const Eigen::Ref<const Eigen::VectorXd>& yv,
    const Eigen::Ref<const Eigen::VectorXd>& eta0,
    const Eigen::Ref<const Eigen::VectorXd>& mu0,
    const Eigen::Ref<const Eigen::VectorXd>& beta_a,
    int rankA,
    const Eigen::Ref<const Eigen::VectorXd>& W0,
    const Eigen::Ref<const Eigen::VectorXd>& z0,
    const Eigen::Ref<const Eigen::MatrixXd>& L_AtWA_in,
    bool L_AtWA_not_ok,
    double dev0,
    int maxit,
    double tol,
    int max_step_halving,
    double tol_scale_rank,
    LogisticCoefAllWorkspaceEigen& ws) {
  const int pX = static_cast<int>(Xm.cols());
  const int p_full = pX + rankA;

  ws.beta.setZero();
  ws.beta.tail(rankA) = beta_a;
  ws.eta = eta0;
  ws.mu = mu0;
  double dev = dev0;

  ws.W = W0;
  ws.z = z0;
  ws.Z.leftCols(pX) = Xm;
  ws.Z.rightCols(rankA) = Am;

  bool converged = false;
  bool has_error = false;
  std::string error_msg;
  int it = 0;
  bool last_used_qr = false;

  const Eigen::Ref<const Eigen::MatrixXd>* L_ptr = L_AtWA_not_ok ? nullptr : &L_AtWA_in;
  auto beta_x_new = ws.beta_new.head(pX);
  auto beta_a_new = ws.beta_new.tail(rankA);

  for (; it < maxit; ++it) {
    SchurStepResult schur_res = wls_schur_step_eigen(
      Xm, Am, ws.z, ws.W,
      ws.sqrtW_schur, ws.Aw_schur, ws.Xw_schur, ws.zw_schur,
      ws.M_schur, ws.B_schur, ws.c_schur, ws.MinvB_schur, ws.Minvc_schur, ws.S_schur, ws.r_schur,
      L_ptr, L_AtWA_not_ok, tol_scale_rank
    );

    if (schur_res.rank_deficient) {
      has_error = true;
      error_msg = "rank_deficient_in_irls";
      break;
    }

    if (schur_res.success) {
      beta_x_new.noalias() = schur_res.beta_x;
      beta_a_new.noalias() = schur_res.beta_a;
      last_used_qr = schur_res.used_qr_for_solve;
    } else {
      int rank_wls = 0;
      bool used_qr = false;
      wls_cholesky_qr_step_eigen(
        ws.Z, ws.z, ws.W,
        ws.sqrtW_wls, ws.Xw_wls, ws.zw_wls,
        ws.XtWX_wls, ws.XtWz_wls,
        ws.beta_new, &rank_wls, &used_qr, tol_scale_rank
      );
      if (rank_wls < p_full) {
        has_error = true;
        error_msg = "rank_deficient_in_irls";
        break;
      }
      last_used_qr = used_qr;
    }

    ws.eta_new.noalias() = Xm * beta_x_new + Am * beta_a_new;
    logistic_mu_stable_eigen(ws.eta_new, ws.mu_new);
    double dev_new = binom_deviance_eigen(yv, ws.eta_new);

    bool has_exact_boundary = (ws.mu_new.array() <= 0.0).any() || (ws.mu_new.array() >= 1.0).any();
    bool need_halving = (!std::isfinite(dev_new)) || (dev_new > dev) || has_exact_boundary;

    if (need_halving) {
      int halving = 0;
      while (halving < max_step_halving) {
        ws.beta_new.noalias() = 0.5 * (ws.beta_new + ws.beta);
        ws.eta_new.noalias() = Xm * beta_x_new + Am * beta_a_new;
        logistic_mu_stable_eigen(ws.eta_new, ws.mu_new);
        dev_new = binom_deviance_eigen(yv, ws.eta_new);
        has_exact_boundary = (ws.mu_new.array() <= 0.0).any() || (ws.mu_new.array() >= 1.0).any();
        if (std::isfinite(dev_new) && dev_new <= dev && !has_exact_boundary) {
          break;
        }
        ++halving;
      }

      if (!std::isfinite(dev_new) || has_exact_boundary) {
        has_error = true;
        error_msg = has_exact_boundary ? "mu_at_boundary_after_halving" : "infinite_deviance_after_halving";
        break;
      }
    }

    double max_diff = (ws.beta_new - ws.beta).cwiseAbs().maxCoeff();
    double rel_dev = std::abs(dev - dev_new) / (0.1 + std::abs(dev));

    ws.beta.noalias() = ws.beta_new;
    ws.eta.noalias() = ws.eta_new;
    ws.mu.noalias() = ws.mu_new;
    dev = dev_new;

    if (rel_dev < tol && max_diff < tol) {
      converged = true;
      break;
    }

    ws.W.noalias() = (ws.mu.array() * (1.0 - ws.mu.array())).matrix();
    ws.z.noalias() = ws.eta + (yv - ws.mu).cwiseQuotient(ws.W);
    L_ptr = nullptr;
    L_AtWA_not_ok = false;
  }

  if (!converged || has_error) {
    ws.beta.setConstant(NA_REAL);
  }

  LogisticCoefAllResult res;
  res.converged = converged;
  res.iter = converged ? (it + 1) : maxit;
  res.deviance = dev;
  res.error = has_error;
  res.error_msg = error_msg;
  res.used_qr_fallback = last_used_qr;
  return res;
}

struct LogisticXStatsWorkspaceEigen {
  Eigen::VectorXd beta_cur;
  Eigen::VectorXd eta_cur;
  Eigen::VectorXd mu_cur;

  Eigen::VectorXd W;
  Eigen::VectorXd z;
  Eigen::VectorXd eta_new;
  Eigen::VectorXd mu_new;
  Eigen::VectorXd beta_new;
  Eigen::MatrixXd Z;

  Eigen::VectorXd sqrtW_schur;
  Eigen::MatrixXd Aw_schur;
  Eigen::MatrixXd Xw_schur;
  Eigen::VectorXd zw_schur;
  Eigen::MatrixXd M_schur;
  Eigen::MatrixXd B_schur;
  Eigen::VectorXd c_schur;
  Eigen::MatrixXd MinvB_schur;
  Eigen::VectorXd Minvc_schur;
  Eigen::MatrixXd S_schur;
  Eigen::VectorXd r_schur;

  Eigen::VectorXd sqrtW_wls;
  Eigen::MatrixXd Xw_wls;
  Eigen::VectorXd zw_wls;
  Eigen::MatrixXd XtWX_wls;
  Eigen::VectorXd XtWz_wls;

  Eigen::MatrixXd Zw_stats;
  Eigen::MatrixXd ZtWZ_stats;
  Eigen::MatrixXd R_Z_stats;
  Eigen::MatrixXd inv_RZ_stats;
  Eigen::VectorXi perm_invZ_stats;
  Eigen::VectorXd se_decomp_stats;
  Eigen::VectorXd se_x_stats;
  Eigen::MatrixXd A_wald;
  Eigen::MatrixXd M_SS_wald;
  Eigen::MatrixXd R_s_iter;
  Eigen::VectorXi perm_s_iter;
  Eigen::VectorXd beta_perm_stats;
  Eigen::VectorXi perm_inv_s_stats;
  Eigen::MatrixXd R_full_iter;
  Eigen::VectorXi perm_full_iter;

  void resize(int n, int rankA, int pX, int p_full, int df_wald) {
    beta_cur.resize(p_full);
    eta_cur.resize(n);
    mu_cur.resize(n);

    W.resize(n);
    z.resize(n);
    eta_new.resize(n);
    mu_new.resize(n);
    beta_new.resize(p_full);
    Z.resize(n, p_full);

    sqrtW_schur.resize(n);
    Aw_schur.resize(n, rankA);
    Xw_schur.resize(n, pX);
    zw_schur.resize(n);
    M_schur.resize(rankA, rankA);
    B_schur.resize(rankA, pX);
    c_schur.resize(rankA);
    MinvB_schur.resize(rankA, pX);
    Minvc_schur.resize(rankA);
    S_schur.resize(pX, pX);
    r_schur.resize(pX);

    sqrtW_wls.resize(n);
    Xw_wls.resize(n, p_full);
    zw_wls.resize(n);
    XtWX_wls.resize(p_full, p_full);
    XtWz_wls.resize(p_full);

    Zw_stats.resize(n, p_full);
    ZtWZ_stats.resize(p_full, p_full);
    R_Z_stats.resize(p_full, p_full);
    inv_RZ_stats.resize(p_full, p_full);
    perm_invZ_stats.resize(p_full);
    se_decomp_stats.resize(p_full);
    se_x_stats.resize(pX);
    A_wald.resize(df_wald, p_full);
    M_SS_wald.resize(df_wald, df_wald);
    R_s_iter.resize(pX, pX);
    perm_s_iter.resize(pX);
    beta_perm_stats.resize(pX);
    perm_inv_s_stats.resize(pX);
    R_full_iter.resize(p_full, p_full);
    perm_full_iter.resize(p_full);
  }
};

static inline double tlstractor_logistic_x_stats(
    const Eigen::Ref<const Eigen::MatrixXd>& Xm,
    int df,
    const Eigen::Ref<const Eigen::MatrixXd>& Am,
    const Eigen::Ref<const Eigen::VectorXd>& yv,
    const Eigen::Ref<const Eigen::VectorXd>& eta0,
    const Eigen::Ref<const Eigen::VectorXd>& mu0,
    const Eigen::Ref<const Eigen::VectorXd>& beta_a0,
    const Eigen::Ref<const Eigen::VectorXd>& W0,
    const Eigen::Ref<const Eigen::VectorXd>& z0,
    const Eigen::Ref<const Eigen::MatrixXd>& L_AtWA_in,
    bool L_AtWA_not_ok,
    double dev0,
    int maxit,
    double tol,
    int max_step_halving,
    double tol_scale_rank,
    bool use_qr_fallback,
    bool refine_W,
    bool* used_qr_fallback_out,
    LogisticXStatsWorkspaceEigen& ws) {
  const int rankA = static_cast<int>(Am.cols());
  const int pX = static_cast<int>(Xm.cols());
  const int p_full = pX + rankA;

  ws.beta_cur.setZero();
  ws.beta_cur.tail(rankA) = beta_a0;
  ws.eta_cur = eta0;
  ws.mu_cur = mu0;
  double dev = dev0;

  ws.W = W0;
  ws.z = z0;
  ws.Z.leftCols(pX) = Xm;
  ws.Z.rightCols(rankA) = Am;

  bool converged = false;
  bool has_error = false;
  std::string error_msg = "";
  int it = 0;
  bool last_used_schur = false;
  bool last_used_wls_qr = false;

  *used_qr_fallback_out = false;

  const Eigen::Ref<const Eigen::MatrixXd>* L_ptr = L_AtWA_not_ok ? nullptr : &L_AtWA_in;

  auto beta_x_new = ws.beta_new.head(pX);
  auto beta_a_new = ws.beta_new.tail(rankA);

  for (; it < maxit; ++it) {
    SchurStepResult schur_res = wls_schur_step_eigen(
      Xm, Am, ws.z, ws.W,
      ws.sqrtW_schur, ws.Aw_schur, ws.Xw_schur, ws.zw_schur,
      ws.M_schur, ws.B_schur, ws.c_schur, ws.MinvB_schur, ws.Minvc_schur, ws.S_schur, ws.r_schur,
      L_ptr, L_AtWA_not_ok, tol_scale_rank, &ws.R_s_iter, &ws.perm_s_iter
    );

    if (schur_res.rank_deficient) {
      has_error = true;
      error_msg = "rank_deficient_in_irls";
      break;
    }

    if (schur_res.success) {
      beta_x_new.noalias() = schur_res.beta_x;
      beta_a_new.noalias() = schur_res.beta_a;
      last_used_schur = true;
      last_used_wls_qr = schur_res.used_qr_for_solve;
    } else {
      int rank_wls = 0;
      bool used_qr = false;
      wls_cholesky_qr_step_eigen(
        ws.Z, ws.z, ws.W,
        ws.sqrtW_wls, ws.Xw_wls, ws.zw_wls,
        ws.XtWX_wls, ws.XtWz_wls,
        ws.beta_new,
        &rank_wls, &used_qr,
        tol_scale_rank, &ws.R_full_iter, &ws.perm_full_iter
      );

      if (rank_wls < static_cast<int>(ws.Z.cols())) {
        has_error = true;
        error_msg = "rank_deficient_in_irls";
        break;
      }
      last_used_schur = false;
      last_used_wls_qr = used_qr;
    }

    ws.eta_new.noalias() = Xm * beta_x_new + Am * beta_a_new;
    logistic_mu_stable_eigen(ws.eta_new, ws.mu_new);
    double dev_new = binom_deviance_eigen(yv, ws.eta_new);

    bool has_exact_boundary = (ws.mu_new.array() <= 0.0).any() || (ws.mu_new.array() >= 1.0).any();
    bool need_halving = (!std::isfinite(dev_new)) || (dev_new > dev) || has_exact_boundary;

    if (need_halving) {
      int halving = 0;
      while (halving < max_step_halving) {
        ws.beta_new.noalias() = 0.5 * (ws.beta_new + ws.beta_cur);
        ws.eta_new.noalias() = Xm * beta_x_new + Am * beta_a_new;
        logistic_mu_stable_eigen(ws.eta_new, ws.mu_new);
        dev_new = binom_deviance_eigen(yv, ws.eta_new);

        has_exact_boundary = (ws.mu_new.array() <= 0.0).any() || (ws.mu_new.array() >= 1.0).any();
        if (std::isfinite(dev_new) && dev_new <= dev && !has_exact_boundary) {
          break;
        }
        ++halving;
      }

      if (!std::isfinite(dev_new) || has_exact_boundary) {
        has_error = true;
        error_msg = has_exact_boundary ? "mu_at_boundary_after_halving" : "infinite_deviance_after_halving";
        break;
      }
    }

    double max_diff = (ws.beta_new - ws.beta_cur).cwiseAbs().maxCoeff();
    double rel_dev = std::abs(dev - dev_new) / (0.1 + std::abs(dev));

    ws.beta_cur.noalias() = ws.beta_new;
    ws.eta_cur.noalias() = ws.eta_new;
    ws.mu_cur.noalias() = ws.mu_new;
    dev = dev_new;

    if (rel_dev < tol && max_diff < tol) {
      converged = true;
      break;
    }

    ws.W.noalias() = (ws.mu_cur.array() * (1.0 - ws.mu_cur.array())).matrix();
    ws.z.noalias() = ws.eta_cur + (yv - ws.mu_cur).cwiseQuotient(ws.W);
    L_ptr = nullptr;
    L_AtWA_not_ok = false;
  }

  auto make_na_stats = [&]() -> double {
    ws.beta_cur.head(pX).setConstant(NA_REAL);
    ws.se_x_stats.setConstant(NA_REAL);
    return NA_REAL;
  };

  auto beta_x_out = ws.beta_cur.head(pX);
  if (!converged) {
    return make_na_stats();
  }

  *used_qr_fallback_out = last_used_wls_qr;

  if (refine_W) {
    ws.W.noalias() = (ws.mu_cur.array() * (1.0 - ws.mu_cur.array())).matrix();
    ws.sqrtW_schur.noalias() = ws.W.array().sqrt().matrix();
    ws.Aw_schur.noalias() = Am;
    ws.Aw_schur.array().colwise() *= ws.sqrtW_schur.array();
    ws.Xw_schur.noalias() = Xm;
    ws.Xw_schur.array().colwise() *= ws.sqrtW_schur.array();

    ws.M_schur.setZero();
    ws.M_schur.template selfadjointView<Eigen::Lower>().rankUpdate(ws.Aw_schur.transpose());

    Eigen::LLT<Eigen::MatrixXd> llt_M(ws.M_schur);
    bool use_cholesky_M = false;
    if (llt_M.info() == Eigen::Success) {
      Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_M.matrixL()).diagonal().cwiseAbs();
      double min_diag = L_diag.minCoeff();
      double max_diag = L_diag.maxCoeff();
      if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
        use_cholesky_M = true;
      }
    }

    if (!use_cholesky_M) {
      if (!use_qr_fallback) {
        *used_qr_fallback_out = true;
        return make_na_stats();
      }

      ws.Zw_stats.leftCols(pX) = ws.Xw_schur;
      ws.Zw_stats.rightCols(rankA) = ws.Aw_schur;

      ws.ZtWZ_stats.setZero();
      ws.ZtWZ_stats.template selfadjointView<Eigen::Lower>().rankUpdate(ws.Zw_stats.transpose());

      bool use_cholesky_Z = false;
      Eigen::LLT<Eigen::MatrixXd> llt_Z(ws.ZtWZ_stats);
      if (llt_Z.info() == Eigen::Success) {
        Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_Z.matrixL()).diagonal().cwiseAbs();
        double min_diag = L_diag.minCoeff();
        double max_diag = L_diag.maxCoeff();
        if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
          ws.R_full_iter = llt_Z.matrixU();
          use_cholesky_Z = true;
        }
      }

      if (!use_cholesky_Z) {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_Z(ws.Zw_stats);
        set_colpiv_qr_threshold(qr_Z, ws.Zw_stats, tol_scale_rank);
        const int rZ = static_cast<int>(qr_Z.rank());
        if (rZ < p_full) {
          *used_qr_fallback_out = true;
          return make_na_stats();
        }
        ws.R_full_iter = qr_Z.matrixR().topLeftCorner(p_full, p_full).template triangularView<Eigen::Upper>();
        ws.perm_full_iter = qr_Z.colsPermutation().indices();
      }

      last_used_schur = false;
      last_used_wls_qr = !use_cholesky_Z;
    } else {
      ws.S_schur.setZero();
      ws.S_schur.template selfadjointView<Eigen::Lower>().rankUpdate(ws.Xw_schur.transpose());
      ws.B_schur.noalias() = ws.Aw_schur.transpose() * ws.Xw_schur;
      ws.MinvB_schur.noalias() = llt_M.solve(ws.B_schur);
      ws.S_schur.noalias() -= ws.B_schur.transpose() * ws.MinvB_schur;

      ws.S_schur = ws.S_schur.template selfadjointView<Eigen::Lower>();
      bool use_cholesky_S = false;

      Eigen::LLT<Eigen::MatrixXd> llt_S(ws.S_schur);
      if (llt_S.info() == Eigen::Success) {
        Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_S.matrixL()).diagonal().cwiseAbs();
        double min_diag = L_diag.minCoeff();
        double max_diag = L_diag.maxCoeff();
        if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
          ws.R_s_iter = llt_S.matrixU();
          use_cholesky_S = true;
        }
      }

      if (!use_cholesky_S) {
        if (!use_qr_fallback) {
          *used_qr_fallback_out = true;
          return make_na_stats();
        }

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_S(ws.S_schur);
        set_colpiv_qr_threshold(qr_S, ws.S_schur, tol_scale_rank);
        const int rS = static_cast<int>(qr_S.rank());
        if (rS < pX) {
          *used_qr_fallback_out = true;
          return make_na_stats();
        }
        ws.R_s_iter = qr_S.matrixR().topLeftCorner(pX, pX).template triangularView<Eigen::Upper>();
        ws.perm_s_iter = qr_S.colsPermutation().indices();
      }

      last_used_schur = true;
      last_used_wls_qr = !use_cholesky_S;
    }
  }

  *used_qr_fallback_out = *used_qr_fallback_out || last_used_wls_qr;

  if (last_used_schur) {
    if (last_used_wls_qr && !use_qr_fallback) {
      return make_na_stats();
    }

    ws.inv_RZ_stats.setIdentity(pX, pX);
    ws.R_s_iter.template triangularView<Eigen::Upper>().solveInPlace(ws.inv_RZ_stats);

    ws.se_x_stats.setZero();
    ws.se_decomp_stats.head(pX) = ws.inv_RZ_stats.topRows(pX).rowwise().squaredNorm().array().sqrt().matrix();
    double wald_stat = NA_REAL;
    if (last_used_wls_qr) {
      for (int j = 0; j < pX; ++j) {
        ws.se_x_stats(ws.perm_s_iter(j)) = ws.se_decomp_stats(j);
      }

      if (df == pX) {
        ws.beta_perm_stats.setZero();
        for (int j = 0; j < pX; ++j) {
          ws.beta_perm_stats(j) = beta_x_out(ws.perm_s_iter(j));
        }
        Eigen::VectorXd R_beta = ws.R_s_iter.triangularView<Eigen::Upper>() * ws.beta_perm_stats;
        wald_stat = R_beta.squaredNorm();
      } else if (df > 0 && df < pX) {
        ws.perm_inv_s_stats.setZero();
        for (int i = 0; i < pX; ++i) {
          ws.perm_inv_s_stats(ws.perm_s_iter(i)) = i;
        }

        ws.A_wald.setZero();
        for (int j = 0; j < df; ++j) {
          int decomp_row = ws.perm_inv_s_stats(j);
          ws.A_wald.row(j) = ws.inv_RZ_stats.row(decomp_row);
        }

        ws.M_SS_wald = (ws.A_wald.topLeftCorner(df, pX) * ws.A_wald.topLeftCorner(df, pX).transpose()).eval();
        Eigen::VectorXd beta_sub = beta_x_out.head(df);
        Eigen::LLT<Eigen::MatrixXd> llt_w(ws.M_SS_wald);
        if (llt_w.info() == Eigen::Success) {
          Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_w.matrixL()).diagonal().cwiseAbs();
          double min_diag = L_diag.minCoeff();
          double max_diag = L_diag.maxCoeff();
          if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
            Eigen::VectorXd tmp = llt_w.solve(beta_sub);
            wald_stat = beta_sub.dot(tmp);
          }
        }
      }
    } else {
      ws.se_x_stats = ws.se_decomp_stats.head(pX);

      if (df == pX) {
        Eigen::VectorXd R_beta = ws.R_s_iter.triangularView<Eigen::Upper>() * beta_x_out;
        wald_stat = R_beta.squaredNorm();
      } else if (df > 0 && df < pX) {
        ws.A_wald.setZero();
        for (int j = 0; j < df; ++j) {
          ws.A_wald.row(j) = ws.inv_RZ_stats.row(j);
        }

        ws.M_SS_wald = (ws.A_wald.topLeftCorner(df, pX) * ws.A_wald.topLeftCorner(df, pX).transpose()).eval();
        Eigen::VectorXd beta_sub = beta_x_out.head(df);
        Eigen::LLT<Eigen::MatrixXd> llt_w(ws.M_SS_wald);
        if (llt_w.info() == Eigen::Success) {
          Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_w.matrixL()).diagonal().cwiseAbs();
          double min_diag = L_diag.minCoeff();
          double max_diag = L_diag.maxCoeff();
          if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
            Eigen::VectorXd tmp = llt_w.solve(beta_sub);
            wald_stat = beta_sub.dot(tmp);
          }
        }
      }
    }

    return wald_stat;
  }

  if (last_used_wls_qr && !use_qr_fallback) {
    return make_na_stats();
  }

  ws.inv_RZ_stats.setIdentity(p_full, p_full);
  ws.R_full_iter.template triangularView<Eigen::Upper>().solveInPlace(ws.inv_RZ_stats);

  ws.se_x_stats.setZero();
  ws.se_decomp_stats.head(p_full) = ws.inv_RZ_stats.topRows(p_full).rowwise().squaredNorm().array().sqrt().matrix();
  double wald_stat = NA_REAL;
  if (last_used_wls_qr) {
    for (int i = 0; i < p_full; ++i) {
      ws.perm_invZ_stats(ws.perm_full_iter(i)) = i;
    }
    for (int j = 0; j < pX; ++j) {
      const int decomp_row = ws.perm_invZ_stats(j);
      ws.se_x_stats(j) = ws.se_decomp_stats(decomp_row);
    }
    if (df > 0 && df <= pX) {
      ws.A_wald.setZero();
      for (int j = 0; j < df; ++j) {
        const int decomp_row = ws.perm_invZ_stats(j);
        ws.A_wald.row(j) = ws.inv_RZ_stats.row(decomp_row);
      }
      ws.M_SS_wald.noalias() = ws.A_wald * ws.A_wald.transpose();
      Eigen::VectorXd beta_sub = beta_x_out.head(df);
      Eigen::LLT<Eigen::MatrixXd> llt_w(ws.M_SS_wald);
      if (llt_w.info() == Eigen::Success) {
        Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_w.matrixL()).diagonal().cwiseAbs();
        double min_diag = L_diag.minCoeff();
        double max_diag = L_diag.maxCoeff();
        if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
          Eigen::VectorXd tmp = llt_w.solve(beta_sub);
          wald_stat = beta_sub.dot(tmp);
        }
      }
    }
  } else {
    ws.se_x_stats = ws.se_decomp_stats.head(pX);
    if (df > 0 && df <= pX) {
      ws.A_wald.setZero();
      for (int j = 0; j < df; ++j) {
        ws.A_wald.row(j) = ws.inv_RZ_stats.row(j);
      }
      ws.M_SS_wald.noalias() = ws.A_wald * ws.A_wald.transpose();
      Eigen::VectorXd beta_sub = beta_x_out.head(df);
      Eigen::LLT<Eigen::MatrixXd> llt_w(ws.M_SS_wald);
      if (llt_w.info() == Eigen::Success) {
        Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_w.matrixL()).diagonal().cwiseAbs();
        double min_diag = L_diag.minCoeff();
        double max_diag = L_diag.maxCoeff();
        if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
          Eigen::VectorXd tmp = llt_w.solve(beta_sub);
          wald_stat = beta_sub.dot(tmp);
        }
      }
    }
  }

  return wald_stat;
}

struct LogisticCoefVcovAWorkspaceEigen {
  Eigen::VectorXd alpha, eta, mu, W, z;
  Eigen::MatrixXd L_AtWA_cached;
  Eigen::VectorXd sqrtW, xw, zw, B, c, MinvB, Minvc;
  Eigen::MatrixXd Aw, M;

  Eigen::MatrixXd Z;
  Eigen::VectorXd z_adj, sqrtW_wls, zw_wls, XtWz_wls, beta_wls;
  Eigen::MatrixXd Xw_wls, XtWX_wls, R_wls_full;
  Eigen::VectorXi perm_wls_full;

  Eigen::VectorXd sqrtW_A, zw_A, XtWz_A;
  Eigen::MatrixXd Xw_A, XtWX_A, R_wls_A;
  Eigen::VectorXi perm_wls_A;

  Eigen::VectorXd alpha_new, eta_new, mu_new;
  Eigen::VectorXd alpha_half;
  Eigen::MatrixXd L_M_final;
  Eigen::VectorXd MinvB_final;
  Eigen::MatrixXd vcov_alpha;
  Eigen::MatrixXd cov_perm;
  Eigen::MatrixXd cov_orig;
  Eigen::MatrixXd cov_full;
  Eigen::MatrixXd Minv;
  bool converged = false;
  int iter = 0;
  double deviance = NA_REAL;
  bool error = false;
  int error_code = 0;
  bool used_qr_fallback = false;

  void resize(int n, int rankA, bool use_qr_fallback) {
    alpha.resize(rankA);
    eta.resize(n);
    mu.resize(n);
    W.resize(n);
    z.resize(n);
    L_AtWA_cached.resize(rankA, rankA);

    sqrtW.resize(n);
    xw.resize(n);
    zw.resize(n);
    B.resize(rankA);
    c.resize(rankA);
    MinvB.resize(rankA);
    Minvc.resize(rankA);
    Aw.resize(n, rankA);
    M.resize(rankA, rankA);

    if (use_qr_fallback) {
      Z.resize(n, rankA + 1);
      z_adj.resize(n);
      sqrtW_wls.resize(n);
      zw_wls.resize(n);
      XtWz_wls.resize(rankA + 1);
      beta_wls.resize(rankA + 1);
      Xw_wls.resize(n, rankA + 1);
      XtWX_wls.resize(rankA + 1, rankA + 1);
      R_wls_full.resize(rankA + 1, rankA + 1);
      perm_wls_full.resize(rankA + 1);

      sqrtW_A.resize(n);
      zw_A.resize(n);
      XtWz_A.resize(rankA);
      Xw_A.resize(n, rankA);
      XtWX_A.resize(rankA, rankA);
      R_wls_A.resize(rankA, rankA);
      perm_wls_A.resize(rankA);
    }

    alpha_new.resize(rankA);
    eta_new.resize(n);
    mu_new.resize(n);
    alpha_half.resize(rankA);
    L_M_final.resize(rankA, rankA);
    MinvB_final.resize(rankA);
    vcov_alpha.resize(rankA, rankA);
    cov_perm.resize(rankA + 1, rankA + 1);
    cov_orig.resize(rankA + 1, rankA + 1);
    cov_full.resize(rankA + 1, rankA + 1);
    Minv.resize(rankA, rankA);
  }
};

static inline bool tlstractor_logistic_coef_vcov_A(
  const Eigen::Ref<const Eigen::VectorXd>& xv,
  const Eigen::Ref<const Eigen::MatrixXd>& Am,
  const Eigen::Ref<const Eigen::VectorXd>& yv,
  const Eigen::Ref<const Eigen::VectorXd>& eta0,
  const Eigen::Ref<const Eigen::VectorXd>& mu0,
  const Eigen::Ref<const Eigen::VectorXd>& beta_a,
  int rankA,
  const Eigen::Ref<const Eigen::MatrixXd>& L_AtWA_in,
  bool L_AtWA_not_ok_precomp,
  const Eigen::Ref<const Eigen::VectorXd>& W0,
  const Eigen::Ref<const Eigen::VectorXd>& z0,
  double dev0,
  LogisticCoefVcovAWorkspaceEigen& ws,
  double beta_x = 0.0,
  bool use_offset = false,
  int maxit = 25,
  double tol = 1e-8,
  int max_step_halving = 10,
  double tol_scale_rank = 1.0,
  bool use_qr_fallback = true,
  bool refine_W = true) {
  constexpr int LOGI_ERR_NONE = 0;
  constexpr int LOGI_ERR_MU_BOUND_INIT = 1;
  constexpr int LOGI_ERR_DEV_INF_INIT = 2;
  constexpr int LOGI_ERR_COLLINEARITY = 3;
  constexpr int LOGI_ERR_RANK_DEFICIENT_IRLS = 4;
  constexpr int LOGI_ERR_MU_BOUND_HALVING = 5;
  constexpr int LOGI_ERR_DEV_INF_HALVING = 6;
  constexpr int LOGI_ERR_CHOLESKY_NO_QR = 7;

  ws.converged = false;
  ws.iter = 0;
  ws.deviance = NA_REAL;
  ws.error = false;
  ws.error_code = LOGI_ERR_NONE;
  ws.used_qr_fallback = false;

  auto fail_with_na = [&](int iter_value, double dev_value, int code) -> bool {
    ws.alpha.setConstant(NA_REAL);
    ws.vcov_alpha.setConstant(NA_REAL);
    ws.converged = false;
    ws.iter = iter_value;
    ws.deviance = dev_value;
    ws.error = true;
    ws.error_code = code;
    return false;
  };

  const int n = static_cast<int>(Am.rows());

  ws.alpha = beta_a;
  Eigen::VectorXd& alpha = ws.alpha;
  
  // Initialize eta, mu, W, z based on whether we can reuse precomp warm start
  Eigen::VectorXd& eta = ws.eta;
  Eigen::VectorXd& mu = ws.mu;
  double dev;
  Eigen::VectorXd& W = ws.W;
  Eigen::VectorXd& z = ws.z;
  
  bool reuse_precomp = !use_offset || (std::abs(beta_x) < LOGISTIC_BETA_X_FLOOR);
  
  if (reuse_precomp) {
    // beta_x effectively 0 and not using fixed offset, reuse precompute state
    eta = eta0;
    mu = mu0;
    dev = dev0;
    W = W0;
    z = z0;
  } else {
    // Initialize with provided beta_x or fixed offset
    eta = eta0 + xv * beta_x;
    mu = Eigen::VectorXd(n);
    logistic_mu_stable_eigen(eta, mu);
    dev = binom_deviance_eigen(yv, eta);

    bool has_exact_boundary = (mu.array() <= 0.0).any() || (mu.array() >= 1.0).any();
    if (has_exact_boundary || !std::isfinite(dev)) {
      return fail_with_na(0, NA_REAL,
                          has_exact_boundary ? LOGI_ERR_MU_BOUND_INIT : LOGI_ERR_DEV_INF_INIT);
    }
    
    // Compute W and z for first iteration
    W.noalias() = (mu.array() * (1.0 - mu.array())).matrix();
    z.noalias() = eta + (yv - mu).cwiseQuotient(W);
  }

  bool converged = false;
  bool has_error = false;
  int error_code = LOGI_ERR_NONE;
  int it = 0;
  bool store_used_qr = false;

  // Extract and set up cached L_AtWA for warm start
  ws.L_AtWA_cached = L_AtWA_in;
  bool L_AtWA_not_ok = L_AtWA_not_ok_precomp && reuse_precomp;
  Eigen::MatrixXd* L_ptr = (reuse_precomp && !L_AtWA_not_ok) ? &ws.L_AtWA_cached : nullptr;
  Eigen::VectorXd& sqrtW = ws.sqrtW;
  Eigen::MatrixXd& Aw = ws.Aw;
  Eigen::VectorXd& xw = ws.xw;
  Eigen::VectorXd& zw = ws.zw;
  Eigen::MatrixXd& M = ws.M;
  Eigen::VectorXd& B = ws.B;
  Eigen::VectorXd& c = ws.c;
  Eigen::VectorXd& MinvB = ws.MinvB;
  Eigen::VectorXd& Minvc = ws.Minvc;

  // WLS fallback work buffers (allocated only if QR fallback is enabled)
  Eigen::MatrixXd& Z = ws.Z;
  Eigen::VectorXd& z_adj = ws.z_adj;
  Eigen::VectorXd& sqrtW_wls = ws.sqrtW_wls;
  Eigen::MatrixXd& Xw_wls = ws.Xw_wls;
  Eigen::VectorXd& zw_wls = ws.zw_wls;
  Eigen::MatrixXd& XtWX_wls = ws.XtWX_wls;
  Eigen::VectorXd& XtWz_wls = ws.XtWz_wls;
  Eigen::VectorXd& beta_wls = ws.beta_wls;
  Eigen::MatrixXd& R_wls_full = ws.R_wls_full;
  Eigen::VectorXi& perm_wls_full = ws.perm_wls_full;

  Eigen::VectorXd& sqrtW_A = ws.sqrtW_A;
  Eigen::MatrixXd& Xw_A = ws.Xw_A;
  Eigen::VectorXd& zw_A = ws.zw_A;
  Eigen::MatrixXd& XtWX_A = ws.XtWX_A;
  Eigen::VectorXd& XtWz_A = ws.XtWz_A;
  Eigen::MatrixXd& R_wls_A = ws.R_wls_A;
  Eigen::VectorXi& perm_wls_A = ws.perm_wls_A;

  if (use_qr_fallback) {
    Z.col(0) = xv;
    Z.rightCols(rankA) = Am;
    perm_wls_full.setLinSpaced(rankA + 1, 0, rankA);
    perm_wls_A.setLinSpaced(rankA, 0, rankA - 1);
  }
  
  Eigen::VectorXd& alpha_new = ws.alpha_new;
  double beta_x_new = beta_x;
  Eigen::VectorXd& eta_new = ws.eta_new;
  Eigen::VectorXd& mu_new = ws.mu_new;
  
  // Variables to save from final iteration for vcov computation
  Eigen::MatrixXd& L_M_final = ws.L_M_final;
  Eigen::VectorXd& MinvB_final = ws.MinvB_final;
  double s_star_final = 0.0;
  bool last_used_qr_fallback = false;
  bool last_fallback_used_qr = false;

  for (; it < maxit; ++it) {
    last_used_qr_fallback = false;
    // Form weighted matrices
    sqrtW.noalias() = W.array().sqrt().matrix();
    Aw.noalias() = Am;
    Aw.array().colwise() *= sqrtW.array();
    xw.noalias() = (xv.array() * sqrtW.array()).matrix();
    zw.noalias() = (z.array() * sqrtW.array()).matrix();
    
    // Cholesky on M
    Eigen::LLT<Eigen::MatrixXd> llt_M;
    bool use_cholesky_M = false;
    
    if (L_ptr != nullptr && it == 0) {
      use_cholesky_M = true;
    } else if (!L_AtWA_not_ok) {
      M.setZero();
      M.template selfadjointView<Eigen::Lower>().rankUpdate(Aw.transpose());
      llt_M.compute(M);
      if (llt_M.info() == Eigen::Success) {
        Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_M.matrixL()).diagonal().cwiseAbs();
        double min_diag = L_diag.minCoeff();
        double max_diag = L_diag.maxCoeff();
        if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
          use_cholesky_M = true;
        }
      }
    }
    
    if (use_cholesky_M) {
      // Compute B = Aw' * xw and c = Aw' * zw
      B.noalias() = Aw.transpose() * xw;
      c.noalias() = Aw.transpose() * zw;
      
      // Define solve functions that use cached L if available
      auto solve_vec = [&](const Eigen::VectorXd& rhs) -> Eigen::VectorXd {
        if (L_ptr) {
          Eigen::VectorXd y = L_ptr->triangularView<Eigen::Lower>().solve(rhs);
          return L_ptr->transpose().triangularView<Eigen::Upper>().solve(y);
        } else {
          return llt_M.solve(rhs);
        }
      };
      
      // Solve M \ B and M \ c
      MinvB = solve_vec(B);
      Minvc = solve_vec(c);
      
      if (use_offset) {
        // Case 2: Fixed beta_x
        // alpha = M^(-1) * (c - B * beta_x)
        alpha_new = Minvc - MinvB * beta_x;
        beta_x_new = beta_x;  // unchanged
      } else {
        // Case 1: Estimate beta_x via Schur complement
        // s_star = x'Wx - B'MinvB (scalar)
        // r = x'Wz - B'Minvc (scalar)
        double s_star = xw.squaredNorm() - B.dot(MinvB);
        
        // Check collinearity
        if (s_star < TOL_COLLINEARITY) {
          return fail_with_na(it, NA_REAL, LOGI_ERR_COLLINEARITY);
        }
      
        double r = xw.dot(zw) - B.dot(Minvc);
        beta_x_new = r / s_star;
        
        // alpha = M^(-1) * (c - B * beta_x)
        alpha_new = Minvc - MinvB * beta_x_new;
        
        // Save for vcov computation
        s_star_final = s_star;
      }
    } else {
      // QR fallback path
      if (use_offset) {
        z_adj.noalias() = z - xv * beta_x;
        int rank_wls = 0;
        bool used_qr = false;
        wls_cholesky_qr_step_eigen(Am, z_adj, W, sqrtW_A, Xw_A, zw_A,
                                   XtWX_A, XtWz_A, alpha_new, &rank_wls, &used_qr,
                                   tol_scale_rank, &R_wls_A, &perm_wls_A);
        if (rank_wls < rankA) {
          return fail_with_na(it, NA_REAL, LOGI_ERR_RANK_DEFICIENT_IRLS);
        }
        beta_x_new = beta_x;
        last_used_qr_fallback = true;
        last_fallback_used_qr = used_qr;
      } else {
        int rank_wls = 0;
        bool used_qr = false;
        wls_cholesky_qr_step_eigen(Z, z, W, sqrtW_wls, Xw_wls, zw_wls,
                                   XtWX_wls, XtWz_wls, beta_wls, &rank_wls, &used_qr,
                                   tol_scale_rank, &R_wls_full, &perm_wls_full);
        if (rank_wls < rankA + 1) {
          return fail_with_na(it, NA_REAL, LOGI_ERR_RANK_DEFICIENT_IRLS);
        }
        beta_x_new = beta_wls(0);
        alpha_new = beta_wls.tail(rankA);
        last_used_qr_fallback = true;
        last_fallback_used_qr = used_qr;
      }
    }
    
    // Save values from final iteration (will be overwritten each iteration)
    if (!last_used_qr_fallback) {
      if (L_ptr != nullptr && it == 0) {
        L_M_final = *L_ptr;
      } else {
        L_M_final = Eigen::MatrixXd(llt_M.matrixL());
      }
      MinvB_final = MinvB;
    }
    
    // Update eta and mu
    eta_new.noalias() = Am * alpha_new + xv * beta_x_new;
    logistic_mu_stable_eigen(eta_new, mu_new);
    double dev_new = binom_deviance_eigen(yv, eta_new);
    
    // Check for boundary or invalid deviance
    bool has_exact_boundary = (mu_new.array() <= 0.0).any() || (mu_new.array() >= 1.0).any();
    bool need_halving = (!std::isfinite(dev_new)) || (dev_new > dev) || has_exact_boundary;
    
    if (need_halving) {
      int halving = 0;
      Eigen::VectorXd& alpha_half = ws.alpha_half;
      alpha_half.noalias() = alpha;
      double beta_x_half = beta_x;
      
      while (halving < max_step_halving) {
        alpha_half = 0.5 * (alpha_new + alpha_half);
        beta_x_half = 0.5 * (beta_x_new + beta_x_half);
        
        eta_new.noalias() = Am * alpha_half + xv * beta_x_half;
        logistic_mu_stable_eigen(eta_new, mu_new);
        dev_new = binom_deviance_eigen(yv, eta_new);
        
        has_exact_boundary = (mu_new.array() <= 0.0).any() || (mu_new.array() >= 1.0).any();
        
        if (std::isfinite(dev_new) && dev_new <= dev && !has_exact_boundary) {
          alpha_new = alpha_half;
          beta_x_new = beta_x_half;
          break;
        }
        ++halving;
      }
      
      if (!std::isfinite(dev_new) || has_exact_boundary) {
        has_error = true;
        error_code = has_exact_boundary ? LOGI_ERR_MU_BOUND_HALVING : LOGI_ERR_DEV_INF_HALVING;
        break;
      }
    }
    
    // Check convergence
    double max_diff_alpha = (alpha_new - alpha).cwiseAbs().maxCoeff();
    double max_diff_beta_x = std::abs(beta_x_new - beta_x);
    double max_diff = std::max(max_diff_alpha, max_diff_beta_x);
    double rel_dev = std::abs(dev - dev_new) / (0.1 + std::abs(dev));
    
    alpha = alpha_new;
    beta_x = beta_x_new;
    eta = eta_new;
    mu = mu_new;
    dev = dev_new;
    
    if (rel_dev < tol && max_diff < tol) {
      converged = true;
      break;
    }

    W.noalias() = (mu.array() * (1.0 - mu.array())).matrix();
    z.noalias() = eta + (yv - mu).cwiseQuotient(W);
    L_ptr = nullptr;
    L_AtWA_not_ok = false;
  }

  store_used_qr = last_fallback_used_qr;
  
  // Optionally recompute W/factorization at final (alpha, beta_x) for vcov stability
  if (converged && !has_error && refine_W) {
    W.noalias() = (mu.array() * (1.0 - mu.array())).matrix();
    z.noalias() = eta + (yv - mu).cwiseQuotient(W);

    sqrtW.noalias() = W.array().sqrt().matrix();
    Aw.noalias() = Am;
    Aw.array().colwise() *= sqrtW.array();

    M.setZero();
    M.template selfadjointView<Eigen::Lower>().rankUpdate(Aw.transpose());

    Eigen::LLT<Eigen::MatrixXd> llt_M_ref(M);
    bool use_cholesky_M_ref = false;
    if (llt_M_ref.info() == Eigen::Success) {
      Eigen::VectorXd L_diag = Eigen::MatrixXd(llt_M_ref.matrixL()).diagonal().cwiseAbs();
      double min_diag = L_diag.minCoeff();
      double max_diag = L_diag.maxCoeff();
      if (min_diag > CHOLESKY_DIAG_FLOOR && max_diag / min_diag < CHOLESKY_COND_LIMIT) {
        use_cholesky_M_ref = true;
      }
    }

    if (use_cholesky_M_ref) {
      L_M_final = Eigen::MatrixXd(llt_M_ref.matrixL());
      last_used_qr_fallback = false;

      if (!use_offset) {
        xw.noalias() = (xv.array() * sqrtW.array()).matrix();
        B.noalias() = Aw.transpose() * xw;
        MinvB = llt_M_ref.solve(B);
        MinvB_final = MinvB;
        s_star_final = xw.squaredNorm() - B.dot(MinvB);
        if (s_star_final < TOL_COLLINEARITY) {
          return fail_with_na(it, NA_REAL, LOGI_ERR_COLLINEARITY);
        }
      }
    } else {
      if (!use_qr_fallback) {
        has_error = true;
        error_code = LOGI_ERR_CHOLESKY_NO_QR;
      } else {
        int rank_wls = 0;
        bool used_qr = false;
        if (use_offset) {
          z_adj.noalias() = z - xv * beta_x;
          wls_cholesky_qr_step_eigen(Am, z_adj, W, sqrtW_A, Xw_A, zw_A,
                                     XtWX_A, XtWz_A, alpha_new, &rank_wls, &used_qr,
                                     tol_scale_rank, &R_wls_A, &perm_wls_A);
          if (rank_wls < rankA) {
            return fail_with_na(it, NA_REAL, LOGI_ERR_RANK_DEFICIENT_IRLS);
          }
          last_used_qr_fallback = true;
          last_fallback_used_qr = used_qr;
        } else {
          wls_cholesky_qr_step_eigen(Z, z, W, sqrtW_wls, Xw_wls, zw_wls,
                                     XtWX_wls, XtWz_wls, beta_wls, &rank_wls, &used_qr,
                                     tol_scale_rank, &R_wls_full, &perm_wls_full);
          if (rank_wls < rankA + 1) {
            return fail_with_na(it, NA_REAL, LOGI_ERR_RANK_DEFICIENT_IRLS);
          }
          last_used_qr_fallback = true;
          last_fallback_used_qr = used_qr;
        }
      }
    }
  }

  store_used_qr = store_used_qr || last_used_qr_fallback;
  
  if (!converged || has_error) {
    ws.alpha.setConstant(NA_REAL);
    ws.vcov_alpha.setConstant(NA_REAL);
    ws.converged = converged;
    ws.iter = converged ? (it + 1) : maxit;
    ws.deviance = dev;
    ws.error = has_error;
    ws.error_code = error_code;
    ws.used_qr_fallback = store_used_qr;
    return false;
  }

  // Compute vcov_alpha
  Eigen::MatrixXd& vcov_alpha = ws.vcov_alpha;

  // vcov_alpha = M^(-1) [+ correction term for estimated beta_x]
  // M^(-1) = L^(-T) * L^(-1)
  if (last_used_qr_fallback) {
    auto build_vcov_from_R = [&](const Eigen::MatrixXd& R,
                                  const Eigen::VectorXi* perm,
                                  int p,
                                  Eigen::Ref<Eigen::MatrixXd> out) {
      out.setZero();
      out.topLeftCorner(p, p).setIdentity();
      Eigen::Ref<Eigen::MatrixXd> cov_perm = out.topLeftCorner(p, p);
      R.transpose().template triangularView<Eigen::Lower>().solveInPlace(cov_perm);
      R.template triangularView<Eigen::Upper>().solveInPlace(cov_perm);
      if (!perm) {
        return;
      }
      ws.cov_orig.topLeftCorner(p, p).setZero();
      for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
          ws.cov_orig((*perm)(i), (*perm)(j)) = cov_perm(i, j);
        }
      }
      cov_perm = ws.cov_orig.topLeftCorner(p, p);
    };

    if (use_offset) {
      build_vcov_from_R(
        R_wls_A,
        last_fallback_used_qr ? &perm_wls_A : nullptr,
        rankA,
        ws.cov_perm
      );
      vcov_alpha = ws.cov_perm.topLeftCorner(rankA, rankA);
    } else {
      build_vcov_from_R(
        R_wls_full,
        last_fallback_used_qr ? &perm_wls_full : nullptr,
        rankA + 1,
        ws.cov_full
      );
      vcov_alpha = ws.cov_full.block(1, 1, rankA, rankA);
    }
  } else {
    Eigen::MatrixXd& Minv = ws.Minv;
    Minv.setIdentity();
    L_M_final.triangularView<Eigen::Lower>().solveInPlace(Minv);
    L_M_final.transpose().triangularView<Eigen::Upper>().solveInPlace(Minv);

    if (use_offset) {
      // Case 2: vcov_alpha = M^(-1)
      vcov_alpha = Minv;
    } else {
      // Case 1: vcov_alpha = M^(-1) + M^(-1)*B*B'*M^(-1) / s_star
      //                    = M^(-1) + MinvB*MinvB' / s_star
      vcov_alpha = Minv + (MinvB_final * MinvB_final.transpose()) / s_star_final;
    }
  }

  // Ensure symmetry
  vcov_alpha = vcov_alpha.template selfadjointView<Eigen::Lower>();

  ws.converged = converged;
  ws.iter = converged ? (it + 1) : maxit;
  ws.deviance = dev;
  ws.error = has_error;
  ws.error_code = LOGI_ERR_NONE;
  ws.used_qr_fallback = store_used_qr;
  return true;
}

// ============================
// Logistic regression for tlstractor
// ============================

// [[Rcpp::export]]
List fill_chunk_with_sumstats_softpass_logistic(List dos_list, List hap_list, IntegerVector idx,
                                                NumericVector sumstats_beta, NumericVector sumstats_se,
                                                List precomp, List control,
                                                bool cond_local, int num_ancs) {
  const int n_snps = idx.size();
  const int n_samples = Rf_nrows(dos_list[0]);
  const int n_laeff = cond_local ? num_ancs - 1 : 0;
  const int n_dos_nonref = num_ancs - 1;
  const int p_x = 1 + n_dos_nonref + n_laeff;

  NumericMatrix beta(n_snps, num_ancs);
  NumericMatrix se(n_snps, num_ancs);
  NumericMatrix zstat(n_snps, num_ancs);
  std::fill(beta.begin(), beta.end(), NA_REAL);
  std::fill(se.begin(), se.end(), NA_REAL);
  std::fill(zstat.begin(), zstat.end(), NA_REAL);
  NumericVector wald(n_snps, NA_REAL);
  LogicalVector used_qr_fallback(n_snps, false);

  NumericMatrix laeff, lase, laz;
  if (cond_local) {
    laeff = NumericMatrix(n_snps, n_laeff);
    lase = NumericMatrix(n_snps, n_laeff);
    laz = NumericMatrix(n_snps, n_laeff);
    std::fill(laeff.begin(), laeff.end(), NA_REAL);
    std::fill(lase.begin(), lase.end(), NA_REAL);
    std::fill(laz.begin(), laz.end(), NA_REAL);
  }

  std::vector<NumericMatrix> dos_mats;
  dos_mats.reserve(n_dos_nonref);
  NumericMatrix z_mat = subset_cols_cast_center_to_numeric(dos_list[0], idx, "dos_list");
  Eigen::Map<Eigen::MatrixXd> Zm(z_mat.begin(), n_samples, n_snps);
  for (int k = 1; k < num_ancs; ++k) {
    dos_mats.emplace_back(subset_cols_cast_center_to_numeric(dos_list[k], idx, "dos_list"));
    Eigen::Map<const Eigen::MatrixXd> dmap(dos_mats.back().begin(), n_samples, n_snps);
    Zm.noalias() += dmap;
  }

  std::vector<NumericMatrix> hap_mats;
  if (cond_local) {
    hap_mats.reserve(n_laeff);
    for (int k = 0; k < n_laeff; ++k) {
      hap_mats.emplace_back(subset_cols_cast_center_to_numeric(hap_list[k], idx, "hap_list"));
    }
  }

  const bool refine_C = control.containsElementNamed("refine_C")
    ? as<bool>(control["refine_C"]) : false;

  NumericMatrix A_r = precomp["A"];
  NumericVector y_r = precomp["y"];
  NumericVector eta_r = precomp["eta"];
  NumericVector mu_r = precomp["mu"];
  NumericVector beta_a_r = precomp["beta"];
  int rankA = as<int>(precomp["rank"]);
  NumericVector W_r = precomp["W"];
  NumericVector z_r = precomp["z"];
  NumericMatrix L_AtWA_r = precomp["L_AtWA"];
  const bool L_AtWA_not_ok = as<bool>(precomp["L_AtWA_not_ok"]);
  const double dev0 = as<double>(precomp["deviance"]);
  NumericMatrix inv_GammaAA_r = precomp["inv_GammaAA"];
  NumericMatrix vcov_A_null_r = precomp["vcov_A_null"];

  Eigen::Map<const Eigen::MatrixXd> Am(A_r.begin(), A_r.nrow(), A_r.ncol());
  Eigen::Map<const Eigen::VectorXd> yv(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> eta0(eta_r.begin(), eta_r.size());
  Eigen::Map<const Eigen::VectorXd> mu0(mu_r.begin(), mu_r.size());
  Eigen::Map<const Eigen::VectorXd> beta_a(beta_a_r.begin(), beta_a_r.size());
  Eigen::Map<const Eigen::VectorXd> W0(W_r.begin(), W_r.size());
  Eigen::Map<const Eigen::VectorXd> z0(z_r.begin(), z_r.size());
  Eigen::Map<const Eigen::MatrixXd> L_AtWA_in(L_AtWA_r.begin(), L_AtWA_r.nrow(), L_AtWA_r.ncol());
  Eigen::Map<const Eigen::MatrixXd> inv_GammaAA_map(inv_GammaAA_r.begin(), inv_GammaAA_r.nrow(), inv_GammaAA_r.ncol());
  Eigen::Map<const Eigen::MatrixXd> vcov_A_null_map(vcov_A_null_r.begin(), vcov_A_null_r.nrow(), vcov_A_null_r.ncol());

  const int p_full = p_x + rankA;
  GMMWorkspace gmm_ws;
  gmm_ws.resize(n_samples, p_full, rankA, false);
  gmm_ws.X_full.rightCols(rankA) = Am;

  LogisticCoefAllWorkspaceEigen coef_ws;
  coef_ws.resize(n_samples, rankA, p_x, p_full);

  Eigen::VectorXd eta_theta(n_samples);
  Eigen::VectorXd mu_theta(n_samples);
  Eigen::VectorXd mu_beta(n_samples);
  Eigen::MatrixXd T_dos = Eigen::MatrixXd::Zero(num_ancs, num_ancs);
  T_dos(0, 0) = 1.0;
  for (int k = 1; k < num_ancs; ++k) {
    T_dos(k, 0) = 1.0;
    T_dos(k, k) = 1.0;
  }

  auto fill_W_matrix = [&](int j) {
    gmm_ws.X_full.col(0).noalias() = Zm.col(j);
    for (int k = 1; k < num_ancs; ++k) {
      const NumericMatrix& dmat = dos_mats[k - 1];
      const double* src_col = dmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
      std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(k).data());
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const NumericMatrix& hmat = hap_mats[k];
        const int col_out = 1 + n_dos_nonref + k;
        const double* src_col = hmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
        std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(col_out).data());
      }
    }
  };

  for (int j = 0; j < n_snps; ++j) {
    bool used_qr_this_snp = false;
    const double curr_beta = sumstats_beta[j];
    const double curr_se = sumstats_se[j];
    const double V_thetaZ = curr_se * curr_se;
    if (!std::isfinite(V_thetaZ)) {
      continue;
    }

    fill_W_matrix(j);
    gmm_ws.Z.noalias() = Zm.col(j);

    LogisticCoefAllResult fit_res = tlstractor_logistic_coef_all(
      gmm_ws.X_full.leftCols(p_x),
      Am,
      yv,
      eta0,
      mu0,
      beta_a,
      rankA,
      W0,
      z0,
      L_AtWA_in,
      L_AtWA_not_ok,
      dev0,
      25,
      1e-8,
      10,
      1.0,
      coef_ws
    );
    if (!fit_res.converged || fit_res.error || coef_ws.beta.size() < p_full || !coef_ws.beta.allFinite()) {
      continue;
    }
    used_qr_this_snp = fit_res.used_qr_fallback;

    gmm_ws.beta_full.noalias() = coef_ws.beta;
    gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;
    logistic_mu_stable_eigen(gmm_ws.X_beta, mu_beta);

    eta_theta.noalias() = eta0;
    eta_theta.noalias() += gmm_ws.Z * curr_beta;
    logistic_mu_stable_eigen(eta_theta, mu_theta);

    if (!delta_opt_gwas_logistic_softpass(yv, gmm_ws.Z, gmm_ws.X_full,
                                  Am,
                                  mu_beta, mu_theta,
                                  mu0,
                                  vcov_A_null_map,
                                  V_thetaZ,
                                  inv_GammaAA_map,
                                  gmm_ws)) {
      continue;
    }
    bool used_qr_tmp = false;
    if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

    direct_fast_build_pseudo_logistic(gmm_ws.Cn, yv, gmm_ws.Z, gmm_ws.X_full,
                                      gmm_ws.X_beta, mu_beta, mu_theta, gmm_ws);
    if (!gmm_ws.ps_XtX.allFinite() || !gmm_ws.ps_Xty.allFinite()) {
      continue;
    }

    used_qr_tmp = false;
    if (!solve_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.ps_Xty, gmm_ws.beta_full, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

    gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;
    logistic_mu_stable_eigen(gmm_ws.X_beta, mu_beta);

    if (refine_C) {
      if (!delta_opt_gwas_logistic_softpass(yv, gmm_ws.Z, gmm_ws.X_full,
                                    Am,
                                    mu_beta, mu_theta,
                                    mu0,
                                    vcov_A_null_map,
                                    V_thetaZ,
                                    inv_GammaAA_map,
                                    gmm_ws)) {
        continue;
      }
      used_qr_tmp = false;
      if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
    }

    direct_fast_build_psXtX_logistic(gmm_ws.Cn, gmm_ws.Z, gmm_ws.X_full,
                     mu_beta, gmm_ws);
    used_qr_tmp = false;
    if (!invert_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.vcov_beta, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
    gmm_ws.vcov_beta /= static_cast<double>(n_samples);
    if (!gmm_ws.vcov_beta.allFinite()) {
      continue;
    }

    Eigen::VectorXd beta_dos_param = gmm_ws.beta_full.head(num_ancs);
    Eigen::MatrixXd cov_dos_param = gmm_ws.vcov_beta.topLeftCorner(num_ancs, num_ancs);
    Eigen::VectorXd beta_dos = T_dos * beta_dos_param;
    Eigen::MatrixXd cov_dos = T_dos * cov_dos_param * T_dos.transpose();
    cov_dos = 0.5 * (cov_dos + cov_dos.transpose());

    for (int k = 0; k < num_ancs; ++k) {
      beta(j, k) = beta_dos(k);
      const double vkk = cov_dos(k, k);
      if (std::isfinite(vkk) && vkk > 0.0) {
        se(j, k) = std::sqrt(vkk);
        const double z = beta(j, k) / se(j, k);
        zstat(j, k) = std::isfinite(z) ? z : NA_REAL;
      }
    }

    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const int col = 1 + n_dos_nonref + k;
        laeff(j, k) = gmm_ws.beta_full(col);
        const double vkk = gmm_ws.vcov_beta(col, col);
        if (std::isfinite(vkk) && vkk > 0.0) {
          lase(j, k) = std::sqrt(vkk);
          const double z = laeff(j, k) / lase(j, k);
          laz(j, k) = std::isfinite(z) ? z : NA_REAL;
        }
      }
    }

    Eigen::VectorXd wald_tmp(num_ancs);
    used_qr_tmp = false;
    if (solve_spd_or_qr_eigen(cov_dos, beta_dos, wald_tmp, &used_qr_tmp)) {
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
      const double wald_stat = beta_dos.dot(wald_tmp);
      if (std::isfinite(wald_stat) && wald_stat >= 0.0) {
        wald[j] = wald_stat;
      }
    }
    used_qr_fallback[j] = used_qr_this_snp;
  }

  List out = List::create(
    Named("beta") = beta,
    Named("se") = se,
    Named("z") = zstat,
    Named("wald") = wald,
    Named("used_qr_fallback") = used_qr_fallback
  );
  if (cond_local) {
    out["laeff"] = laeff;
    out["lase"] = lase;
    out["laz"] = laz;
  }
  return out;
}

// [[Rcpp::export]]
List fill_chunk_with_sumstats_softfail_logistic(List dos_list, List hap_list, IntegerVector idx,
                                                NumericVector sumstats_beta, NumericVector sumstats_se,
                                                List precomp, List control,
                                                bool cond_local, int num_ancs) {
  const int n_snps = idx.size();
  const int n_samples = Rf_nrows(dos_list[0]);
  const int n_laeff = cond_local ? num_ancs - 1 : 0;
  const int n_dos_nonref = num_ancs - 1;
  const int p_x = 1 + n_dos_nonref + n_laeff;

  NumericMatrix beta(n_snps, num_ancs);
  NumericMatrix se(n_snps, num_ancs);
  NumericMatrix zstat(n_snps, num_ancs);
  std::fill(beta.begin(), beta.end(), NA_REAL);
  std::fill(se.begin(), se.end(), NA_REAL);
  std::fill(zstat.begin(), zstat.end(), NA_REAL);
  NumericVector wald(n_snps, NA_REAL);
  LogicalVector used_qr_fallback(n_snps, false);

  NumericMatrix laeff, lase, laz;
  if (cond_local) {
    laeff = NumericMatrix(n_snps, n_laeff);
    lase = NumericMatrix(n_snps, n_laeff);
    laz = NumericMatrix(n_snps, n_laeff);
    std::fill(laeff.begin(), laeff.end(), NA_REAL);
    std::fill(lase.begin(), lase.end(), NA_REAL);
    std::fill(laz.begin(), laz.end(), NA_REAL);
  }

  std::vector<NumericMatrix> dos_mats;
  dos_mats.reserve(n_dos_nonref);
  NumericMatrix z_mat = subset_cols_cast_center_to_numeric(dos_list[0], idx, "dos_list");
  Eigen::Map<Eigen::MatrixXd> Zm(z_mat.begin(), n_samples, n_snps);
  for (int k = 1; k < num_ancs; ++k) {
    dos_mats.emplace_back(subset_cols_cast_center_to_numeric(dos_list[k], idx, "dos_list"));
    Eigen::Map<const Eigen::MatrixXd> dmap(dos_mats.back().begin(), n_samples, n_snps);
    Zm.noalias() += dmap;
  }

  std::vector<NumericMatrix> hap_mats;
  if (cond_local) {
    hap_mats.reserve(n_laeff);
    for (int k = 0; k < n_laeff; ++k) {
      hap_mats.emplace_back(subset_cols_cast_center_to_numeric(hap_list[k], idx, "hap_list"));
    }
  }

  const bool refine_C = control.containsElementNamed("refine_C")
    ? as<bool>(control["refine_C"]) : false;
  const bool use_offset = control.containsElementNamed("use_offset")
    ? as<bool>(control["use_offset"]) : false;
  const bool use_qr_fallback = control.containsElementNamed("use_qr_fallback")
    ? as<bool>(control["use_qr_fallback"]) : true;
  const bool refine_W = control.containsElementNamed("refine_W")
    ? as<bool>(control["refine_W"]) : false;

  NumericMatrix A_r = precomp["A"];
  NumericVector y_r = precomp["y"];
  NumericVector eta_r = precomp["eta"];
  NumericVector mu_r = precomp["mu"];
  NumericVector beta_a_r = precomp["beta"];
  int rankA = as<int>(precomp["rank"]);
  NumericVector W_r = precomp["W"];
  NumericVector z_r = precomp["z"];
  NumericMatrix L_AtWA_r = precomp["L_AtWA"];
  const bool L_AtWA_not_ok = as<bool>(precomp["L_AtWA_not_ok"]);
  const double dev0 = as<double>(precomp["deviance"]);

  Eigen::Map<const Eigen::MatrixXd> Am(A_r.begin(), A_r.nrow(), A_r.ncol());
  Eigen::Map<const Eigen::VectorXd> yv(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> eta0(eta_r.begin(), eta_r.size());
  Eigen::Map<const Eigen::VectorXd> mu0(mu_r.begin(), mu_r.size());
  Eigen::Map<const Eigen::VectorXd> beta_a(beta_a_r.begin(), beta_a_r.size());
  Eigen::Map<const Eigen::VectorXd> W0(W_r.begin(), W_r.size());
  Eigen::Map<const Eigen::VectorXd> z0(z_r.begin(), z_r.size());
  Eigen::Map<const Eigen::MatrixXd> L_AtWA_in(L_AtWA_r.begin(), L_AtWA_r.nrow(), L_AtWA_r.ncol());

  const int p_full = p_x + rankA;
  GMMWorkspace gmm_ws;
  gmm_ws.resize(n_samples, p_full, rankA, true);
  gmm_ws.X_full.rightCols(rankA) = Am;

  LogisticCoefAllWorkspaceEigen coef_ws;
  coef_ws.resize(n_samples, rankA, p_x, p_full);

  LogisticCoefVcovAWorkspaceEigen coef_A_ws;
  coef_A_ws.resize(n_samples, rankA, use_qr_fallback);

  Eigen::VectorXd eta_theta(n_samples);
  Eigen::VectorXd mu_theta(n_samples);
  Eigen::VectorXd mu_beta(n_samples);
  Eigen::MatrixXd T_dos = Eigen::MatrixXd::Zero(num_ancs, num_ancs);
  T_dos(0, 0) = 1.0;
  for (int k = 1; k < num_ancs; ++k) {
    T_dos(k, 0) = 1.0;
    T_dos(k, k) = 1.0;
  }

  auto fill_W_matrix = [&](int j) {
    gmm_ws.X_full.col(0).noalias() = Zm.col(j);
    for (int k = 1; k < num_ancs; ++k) {
      const NumericMatrix& dmat = dos_mats[k - 1];
      const double* src_col = dmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
      std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(k).data());
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const NumericMatrix& hmat = hap_mats[k];
        const int col_out = 1 + n_dos_nonref + k;
        const double* src_col = hmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
        std::copy(src_col, src_col + n_samples, gmm_ws.X_full.col(col_out).data());
      }
    }
  };

  for (int j = 0; j < n_snps; ++j) {
    bool used_qr_this_snp = false;
    const double curr_beta = sumstats_beta[j];
    const double curr_se = sumstats_se[j];
    const double V_thetaZ = curr_se * curr_se;
    if (!std::isfinite(V_thetaZ)) {
      continue;
    }

    fill_W_matrix(j);
    gmm_ws.Z.noalias() = Zm.col(j);

    LogisticCoefAllResult fit_res = tlstractor_logistic_coef_all(
      gmm_ws.X_full.leftCols(p_x),
      Am,
      yv,
      eta0,
      mu0,
      beta_a,
      rankA,
      W0,
      z0,
      L_AtWA_in,
      L_AtWA_not_ok,
      dev0,
      25,
      1e-8,
      10,
      1.0,
      coef_ws
    );
    if (!fit_res.converged || fit_res.error || coef_ws.beta.size() < p_full || !coef_ws.beta.allFinite()) {
      continue;
    }
    used_qr_this_snp = fit_res.used_qr_fallback;

    gmm_ws.beta_full.noalias() = coef_ws.beta;
    gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;
    logistic_mu_stable_eigen(gmm_ws.X_beta, mu_beta);

    bool delta_ok = false;
    bool coefA_ok = tlstractor_logistic_coef_vcov_A(
      gmm_ws.Z,
      Am,
      yv,
      eta0,
      mu0,
      beta_a,
      rankA,
      L_AtWA_in,
      L_AtWA_not_ok,
      W0,
      z0,
      dev0,
      coef_A_ws,
      curr_beta,
      use_offset,
      25,
      1e-8,
      10,
      1.0,
      use_qr_fallback,
      refine_W
    );
    if (!coefA_ok || coef_A_ws.error || !coef_A_ws.converged || !coef_A_ws.vcov_alpha.allFinite()) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || coef_A_ws.used_qr_fallback;

    eta_theta.noalias() = Am * coef_A_ws.alpha + gmm_ws.Z * curr_beta;
    logistic_mu_stable_eigen(eta_theta, mu_theta);

    bool used_qr_delta = false;
    delta_ok = delta_opt_gwas_logistic_softfail(
      yv,
      gmm_ws.Z,
      gmm_ws.X_full,
      Am,
      mu_beta,
      mu_theta,
      coef_A_ws.vcov_alpha,
      V_thetaZ,
      use_offset,
      gmm_ws,
      &used_qr_delta
    );
    if (!delta_ok) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_delta;

    bool used_qr_tmp = false;
    if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

    direct_fast_build_pseudo_logistic(gmm_ws.Cn, yv, gmm_ws.Z, gmm_ws.X_full,
                      gmm_ws.X_beta, mu_beta, mu_theta, gmm_ws);
    if (!gmm_ws.ps_XtX.allFinite() || !gmm_ws.ps_Xty.allFinite()) {
      continue;
    }

    used_qr_tmp = false;
    if (!solve_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.ps_Xty, gmm_ws.beta_full, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;

    gmm_ws.X_beta.noalias() = gmm_ws.X_full * gmm_ws.beta_full;
    logistic_mu_stable_eigen(gmm_ws.X_beta, mu_beta);

    if (refine_C) {
      used_qr_delta = false;
      if (!delta_opt_gwas_logistic_softfail(
            yv,
            gmm_ws.Z,
            gmm_ws.X_full,
            Am,
            mu_beta,
            mu_theta,
            coef_A_ws.vcov_alpha,
            V_thetaZ,
            use_offset,
        gmm_ws,
        &used_qr_delta)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_delta;
      used_qr_tmp = false;
      if (!invert_spd_or_qr_eigen(gmm_ws.inv_C, gmm_ws.Cn, &used_qr_tmp)) {
        continue;
      }
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
    }

    direct_fast_build_psXtX_logistic(gmm_ws.Cn, gmm_ws.Z, gmm_ws.X_full,
                                     mu_beta, gmm_ws);
    used_qr_tmp = false;
    if (!invert_spd_or_qr_eigen(gmm_ws.ps_XtX, gmm_ws.vcov_beta, &used_qr_tmp)) {
      continue;
    }
    used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
    gmm_ws.vcov_beta /= static_cast<double>(n_samples);
    if (!gmm_ws.vcov_beta.allFinite()) {
      continue;
    }

    Eigen::VectorXd beta_dos_param = gmm_ws.beta_full.head(num_ancs);
    Eigen::MatrixXd cov_dos_param = gmm_ws.vcov_beta.topLeftCorner(num_ancs, num_ancs);
    Eigen::VectorXd beta_dos = T_dos * beta_dos_param;
    Eigen::MatrixXd cov_dos = T_dos * cov_dos_param * T_dos.transpose();
    cov_dos = 0.5 * (cov_dos + cov_dos.transpose());

    for (int k = 0; k < num_ancs; ++k) {
      beta(j, k) = beta_dos(k);
      const double vkk = cov_dos(k, k);
      if (std::isfinite(vkk) && vkk > 0.0) {
        se(j, k) = std::sqrt(vkk);
        const double z = beta(j, k) / se(j, k);
        zstat(j, k) = std::isfinite(z) ? z : NA_REAL;
      }
    }

    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const int col = 1 + n_dos_nonref + k;
        laeff(j, k) = gmm_ws.beta_full(col);
        const double vkk = gmm_ws.vcov_beta(col, col);
        if (std::isfinite(vkk) && vkk > 0.0) {
          lase(j, k) = std::sqrt(vkk);
          const double z = laeff(j, k) / lase(j, k);
          laz(j, k) = std::isfinite(z) ? z : NA_REAL;
        }
      }
    }

    Eigen::VectorXd wald_tmp(num_ancs);
    used_qr_tmp = false;
    if (solve_spd_or_qr_eigen(cov_dos, beta_dos, wald_tmp, &used_qr_tmp)) {
      used_qr_this_snp = used_qr_this_snp || used_qr_tmp;
      const double wald_stat = beta_dos.dot(wald_tmp);
      if (std::isfinite(wald_stat) && wald_stat >= 0.0) {
        wald[j] = wald_stat;
      }
    }
    used_qr_fallback[j] = used_qr_this_snp;
  }

  List out = List::create(
    Named("beta") = beta,
    Named("se") = se,
    Named("z") = zstat,
    Named("wald") = wald,
    Named("used_qr_fallback") = used_qr_fallback
  );
  if (cond_local) {
    out["laeff"] = laeff;
    out["lase"] = lase;
    out["laz"] = laz;
  }
  return out;
}

// [[Rcpp::export]]
List fill_chunk_without_sumstats_logistic(List dos_list, List hap_list, IntegerVector idx,
                                          List precomp, List control, bool cond_local, int num_ancs) {
  const int n_snps = idx.size();
  const int n_samples = Rf_nrows(dos_list[0]);
  const int n_laeff = cond_local ? num_ancs - 1 : 0;
  const int p_total = num_ancs + n_laeff;

  NumericMatrix beta(n_snps, num_ancs);
  NumericMatrix se(n_snps, num_ancs);
  NumericMatrix zstat(n_snps, num_ancs);
  std::fill(beta.begin(), beta.end(), NA_REAL);
  std::fill(se.begin(), se.end(), NA_REAL);
  std::fill(zstat.begin(), zstat.end(), NA_REAL);
  NumericVector wald(n_snps, NA_REAL);
  LogicalVector used_qr_fallback(n_snps, false);

  NumericMatrix laeff, lase, laz;
  if (cond_local) {
    laeff = NumericMatrix(n_snps, n_laeff);
    lase = NumericMatrix(n_snps, n_laeff);
    laz = NumericMatrix(n_snps, n_laeff);
    std::fill(laeff.begin(), laeff.end(), NA_REAL);
    std::fill(lase.begin(), lase.end(), NA_REAL);
    std::fill(laz.begin(), laz.end(), NA_REAL);
  }

  std::vector<NumericMatrix> dos_mats;
  dos_mats.reserve(num_ancs);
  for (int k = 0; k < num_ancs; ++k) {
    dos_mats.emplace_back(subset_cols_cast_center_to_numeric(dos_list[k], idx, "dos_list"));
  }

  std::vector<NumericMatrix> hap_mats;
  if (cond_local) {
    hap_mats.reserve(num_ancs);
    for (int k = 0; k < num_ancs; ++k) {
      hap_mats.emplace_back(subset_cols_cast_center_to_numeric(hap_list[k], idx, "hap_list"));
    }
  }

  const bool use_qr_fallback = control.containsElementNamed("use_qr_fallback")
    ? as<bool>(control["use_qr_fallback"]) : true;
  const bool refine_W = control.containsElementNamed("refine_W")
    ? as<bool>(control["refine_W"]) : false;

  NumericMatrix A_r = precomp["A"];
  NumericVector y_r = precomp["y"];
  NumericVector eta_r = precomp["eta"];
  NumericVector mu_r = precomp["mu"];
  NumericVector beta_a_r = precomp["beta"];
  NumericVector W_r = precomp["W"];
  NumericVector z_r = precomp["z"];
  NumericMatrix L_AtWA_r = precomp["L_AtWA"];
  const bool L_AtWA_not_ok = as<bool>(precomp["L_AtWA_not_ok"]);
  const double dev0 = as<double>(precomp["deviance"]);

  Eigen::Map<const Eigen::MatrixXd> Am(A_r.begin(), A_r.nrow(), A_r.ncol());
  Eigen::Map<const Eigen::VectorXd> yv(y_r.begin(), y_r.size());
  Eigen::Map<const Eigen::VectorXd> eta0(eta_r.begin(), eta_r.size());
  Eigen::Map<const Eigen::VectorXd> mu0(mu_r.begin(), mu_r.size());
  Eigen::Map<const Eigen::VectorXd> beta_a0(beta_a_r.begin(), beta_a_r.size());
  Eigen::Map<const Eigen::VectorXd> W0(W_r.begin(), W_r.size());
  Eigen::Map<const Eigen::VectorXd> z0(z_r.begin(), z_r.size());
  Eigen::Map<const Eigen::MatrixXd> L_AtWA_in(L_AtWA_r.begin(), L_AtWA_r.nrow(), L_AtWA_r.ncol());

  const int rankA = as<int>(precomp["rank"]);
  const int pX = p_total;
  const int p_full = pX + rankA;
  const int df_wald = num_ancs;
  LogisticXStatsWorkspaceEigen logistic_workspace;
  logistic_workspace.resize(n_samples, rankA, pX, p_full, df_wald);

  NumericMatrix X(n_samples, p_total);
  Eigen::Map<Eigen::MatrixXd> Xm(X.begin(), n_samples, p_total);
  for (int j = 0; j < n_snps; ++j) {
    for (int k = 0; k < num_ancs; ++k) {
      const NumericMatrix& dmat = dos_mats[k];
      const double* src_col = dmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
      double* dst_col = X.begin() + static_cast<R_xlen_t>(k) * n_samples;
      std::copy(src_col, src_col + n_samples, dst_col);
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const NumericMatrix& hmat = hap_mats[k];
        const int col_out = num_ancs + k;
        const double* src_col = hmat.begin() + static_cast<R_xlen_t>(j) * n_samples;
        double* dst_col = X.begin() + static_cast<R_xlen_t>(col_out) * n_samples;
        std::copy(src_col, src_col + n_samples, dst_col);
      }
    }

    bool used_qr_flag = false;
    const double wj = tlstractor_logistic_x_stats(
      Xm,
      num_ancs,
      Am,
      yv,
      eta0,
      mu0,
      beta_a0,
      W0,
      z0,
      L_AtWA_in,
      L_AtWA_not_ok,
      dev0,
      25,
      1e-8,
      10,
      1.0,
      use_qr_fallback,
      refine_W,
      &used_qr_flag,
      logistic_workspace
    );
    used_qr_fallback[j] = used_qr_flag;

    auto b = logistic_workspace.beta_cur.head(pX);
    const Eigen::VectorXd& s = logistic_workspace.se_x_stats;
    for (int k = 0; k < num_ancs; ++k) {
      beta(j, k) = b(k);
      se(j, k) = s(k);
    }
    if (cond_local) {
      for (int k = 0; k < n_laeff; ++k) {
        const int col_in = num_ancs + k;
        laeff(j, k) = b(col_in);
        lase(j, k) = s(col_in);
      }
    }

    if (std::isfinite(wj) && wj >= 0.0) {
      wald[j] = wj;
    }
  }

  Eigen::Map<const Eigen::ArrayXXd> beta_arr(beta.begin(), n_snps, num_ancs);
  Eigen::Map<const Eigen::ArrayXXd> se_arr(se.begin(), n_snps, num_ancs);
  Eigen::Map<Eigen::ArrayXXd> z_arr(zstat.begin(), n_snps, num_ancs);
  z_arr = beta_arr / se_arr;
  if (cond_local) {
    Eigen::Map<const Eigen::ArrayXXd> laeff_arr(laeff.begin(), n_snps, n_laeff);
    Eigen::Map<const Eigen::ArrayXXd> lase_arr(lase.begin(), n_snps, n_laeff);
    Eigen::Map<Eigen::ArrayXXd> laz_arr(laz.begin(), n_snps, n_laeff);
    laz_arr = laeff_arr / lase_arr;
  }

  List out = List::create(
    Named("beta") = beta,
    Named("se") = se,
    Named("z") = zstat,
    Named("wald") = wald,
    Named("used_qr_fallback") = used_qr_fallback
  );
  if (cond_local) {
    out["laeff"] = laeff;
    out["lase"] = lase;
    out["laz"] = laz;
  }
  return out;
}