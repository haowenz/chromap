#ifndef MAPPING_H_
#define MAPPING_H_

namespace chromap {

// An interface for various mapping formats.
class Mapping {
 public:
  Mapping() {}
  //// Defines the orders of mappings. Sort by mapping positions first, then
  //// sorted by barcode and other fields if available. Make sure to consider
  //// enough field so that the order is always deterministic.
  // virtual bool operator<(const Mapping &m) const = 0;
  //// Return true if two mappings are the same.
  // virtual bool operator==(const Mapping &m) const = 0;
  // Return true if two mappings are the same. For paired-end mappings, return
  // true if the mapping intervals are the same. For single-end mappings, return
  // true if the 5' mapping positions and strands are the same. This is
  // different from the previous function as this function does not require
  // barcodes to be the same.
  // virtual bool IsSamePosition(const Mapping &m) const = 0;
  // Return true if the mapping strand is positive. For paired-end reads, return
  // true if the mapping strand of the first read is positive.
  virtual bool IsPositiveStrand() const = 0;
  // Barcodes are encoded by 64-bit integers. This function will return the
  // encoded barcode. For mapping without barcode, this function will return 0;
  virtual uint64_t GetBarcode() const = 0;
  // Return inclusive mapping start position.
  virtual uint32_t GetStartPosition() const = 0;
  // Return exclusive mapping start position.
  virtual uint32_t GetEndPosition() const = 0;
  //// Return the total byte size of the mapping data field.
  // virtual uint16_t GetByteSize() const = 0;
  //// Write this mapping to a temp mapping output file in binary.
  // virtual size_t WriteToFile(FILE *temp_mapping_output_file) const = 0;
  //// Load this mapping fomr a temp mapping output file.
  // virtual size_t LoadFromFile(FILE *temp_mapping_output_file) = 0;
  // Perform Tn5 shift, which will change the start and end positions. Note that
  // currently this can only be done on the mappings that are represented by
  // intervals.
  virtual void Tn5Shift() = 0;
};

}  // namespace chromap

#endif  // MAPPING_H_
