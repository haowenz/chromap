#ifndef SUMMARY_METADATA_H_
#define SUMMARY_METADATA_H_

#include <string>

#include <stdio.h>
#include <stdint.h>
#include <cmath>

#include "khash.h"
#include "utils.h"

// The class summarizes the overall mapping metadata 

namespace chromap {

enum SummaryMetadataField {
  SUMMARY_METADATA_TOTAL = 0,
  SUMMARY_METADATA_DUP,
  SUMMARY_METADATA_MAPPED,
  SUMMARY_METADATA_LOWMAPQ,
	SUMMARY_METADATA_CACHEHIT,
  SUMMARY_METADATA_CARDINALITY,
  SUMMARY_METADATA_FIELDS
};

struct _barcodeSummaryMetadata {
  int counts[SUMMARY_METADATA_FIELDS];
  _barcodeSummaryMetadata() {
    memset(counts, 0, sizeof(int) * SUMMARY_METADATA_FIELDS);
  }
};


KHASH_MAP_INIT_INT64(k64_barcode_metadata, struct _barcodeSummaryMetadata)

class SummaryMetadata {
 public:
  SummaryMetadata() {
    barcode_metadata_ = kh_init(k64_barcode_metadata);
    barcode_length_ = 16;
  }
  ~SummaryMetadata() {
    kh_destroy(k64_barcode_metadata, barcode_metadata_);
  }

  inline double inverse_logit(double frip) {
    return (1.0/(1.0 + std::exp(-frip)));
  }

  inline void OutputCounts(const char *barcode, const int *counts, FILE *fp, std::vector<double> frip_est_coeffs, bool output_num_cache_slots_info)
  {
    // define variables to store values
    size_t num_total = counts[SUMMARY_METADATA_TOTAL];
    size_t num_dup = counts[SUMMARY_METADATA_DUP]; 
    
    size_t num_mapped = counts[SUMMARY_METADATA_MAPPED];
    size_t num_unmapped = num_total - num_mapped;

    size_t num_lowmapq = counts[SUMMARY_METADATA_LOWMAPQ];
    size_t num_cachehit = counts[SUMMARY_METADATA_CACHEHIT];
    double fric = (num_mapped != 0) ? (double) num_cachehit / (double) num_mapped : 0.0;

    size_t num_cache_slots = counts[SUMMARY_METADATA_CARDINALITY];

    // compute the estimated frip
    double est_frip = (fric != 0.0) ? inverse_logit(frip_est_coeffs[0] + /* constant */
                                           (frip_est_coeffs[1] * fric) +
                                           (frip_est_coeffs[2] * num_dup) +
                                           (frip_est_coeffs[3] * num_unmapped)  +
                                           (frip_est_coeffs[4] * num_lowmapq)) : 0.0;

    // print out data for current barcode
    if (!output_num_cache_slots_info) {
      fprintf(fp, "%s,%ld,%ld,%ld,%ld,%ld,%.5lf,%.5lf\n", 
              barcode,
              num_total,
              num_dup,
              num_unmapped,
              num_lowmapq,
              num_cachehit,
              fric,
              est_frip);
    } else {
      fprintf(fp, "%s,%ld,%ld,%ld,%ld,%ld,%.5lf,%.5lf,%ld\n", 
              barcode,
              num_total,
              num_dup,
              num_unmapped,
              num_lowmapq,
              num_cachehit,
              fric,
              est_frip,
              num_cache_slots);
    }
  }

  void Output(const char *filename, bool has_white_list, std::vector<double> frip_est_coeffs, bool output_num_cache_slots_info) {
    FILE *fp = fopen(filename, "w");

    // Change summary file header depending on options
    if (!output_num_cache_slots_info)
      fprintf(fp, "barcode,total,duplicate,unmapped,lowmapq,cachehit,fric,estfrip\n");   
    else
      fprintf(fp, "barcode,total,duplicate,unmapped,lowmapq,cachehit,fric,estfrip,numcacheslots\n"); 

    khiter_t k;
    for (k = kh_begin(barcode_metadata_); k != kh_end(barcode_metadata_); ++k)
      if (kh_exist(barcode_metadata_, k)) {
        OutputCounts(
                    Seed2Sequence(kh_key(barcode_metadata_, k), barcode_length_).c_str(),
                    kh_value(barcode_metadata_, k).counts, 
                    fp,
                    frip_est_coeffs,
                    output_num_cache_slots_info
                    );
      }
    if (has_white_list) {
      OutputCounts(
                   "non-whitelist", 
                   nonwhitelist_summary_.counts, 
                   fp,
                   frip_est_coeffs,
                   output_num_cache_slots_info
                   ) ;
    }
    fclose(fp);
  }

  void UpdateCount(uint64_t barcode, int type, int change) {
    int khash_return_code;
    khiter_t barcode_metadata_iter = kh_put(k64_barcode_metadata, barcode_metadata_, barcode, &khash_return_code);
    if (khash_return_code) {
      struct _barcodeSummaryMetadata nb;
      kh_value(barcode_metadata_, barcode_metadata_iter) = nb;
    }
    kh_value(barcode_metadata_, barcode_metadata_iter).counts[type] += change;
  }

  void UpdateNonWhitelistCount(int type, int change) {
    nonwhitelist_summary_.counts[type] += change;
  }

  void SetBarcodeLength(int l) {
    barcode_length_ = l;
  }

  // In SAM format for paired-end data, some count will be counted twice
  void AdjustPairedEndOverCount() {
    khiter_t k;
    for (k = kh_begin(barcode_metadata_); k != kh_end(barcode_metadata_); ++k)
      if (kh_exist(barcode_metadata_, k)) {
        kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_DUP] /= 2 ;
        kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_LOWMAPQ] /= 2 ;
        kh_value(barcode_metadata_, k).counts[SUMMARY_METADATA_MAPPED] /= 2 ;
      } 
  }

 private:
  khash_t(k64_barcode_metadata) *barcode_metadata_;    
  struct _barcodeSummaryMetadata nonwhitelist_summary_;  // summarize the fragments with no barcode information 
  int barcode_length_;

  std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask_ = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(
          Uint8ToChar((seed >> ((seed_length - 1 - i) * 2)) & mask_));
    }
    return sequence;
  }
};

} // namespace chromap

#endif
