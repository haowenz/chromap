#ifndef TEMPMAPPING_H_
#define TEMPMAPPING_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "bed_mapping.h"
#include "mapping.h"
#include "paf_mapping.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"

namespace chromap {

template <typename MappingRecord>
struct TempMappingFileHandle {
  std::string file_path;
  FILE* file;
  uint32_t num_mappings;
  uint32_t block_size;
  uint32_t current_rid;
  uint32_t current_mapping_index;
  uint32_t num_mappings_on_current_rid;
  uint32_t num_loaded_mappings_on_current_rid;
  // This vector only keep mappings on the same ref seq.
  std::vector<MappingRecord> mappings;

  inline const MappingRecord& GetCurrentMapping() const {
    return mappings[current_mapping_index];
  }

  inline bool HasMappings() const { return num_mappings != 0; }

  inline void InitializeTempMappingLoading(uint32_t temp_mapping_block_size) {
    file = fopen(file_path.c_str(), "rb");
    if (file == NULL) {
      std::cerr << "Temporary file " << file_path << " is missing.\n" ;
    }
    assert(file != NULL);
    num_mappings = 0;
    block_size = temp_mapping_block_size;
    current_rid = 0;
    current_mapping_index = 0;
    fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
    num_loaded_mappings_on_current_rid = 0;
    mappings.resize(block_size);
    // std::cerr << "Block size: " << block_size << ", initialize temp file " <<
    // file_path << "\n";
  }

  inline void FinalizeTempMappingLoading() { fclose(file); }

  inline void LoadTempMappingBlock(uint32_t num_reference_sequences) {
    num_mappings = 0;
    while (num_mappings == 0) {
      // Only keep mappings on one ref seq, which means # mappings in buffer can
      // be less than block size Two cases: current ref seq has remainings or
      // not
      if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
        // Check if # remains larger than block size
        uint32_t num_mappings_to_load_on_current_rid =
            num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
        if (num_mappings_to_load_on_current_rid > block_size) {
          num_mappings_to_load_on_current_rid = block_size;
        }
        // std::cerr << num_mappings_to_load_on_current_rid << " " <<
        // num_loaded_mappings_on_current_rid << " " <<
        // num_mappings_on_current_rid << "\n"; std::cerr << mappings.size() <<
        // "\n";
        fread(mappings.data(), sizeof(MappingRecord),
              num_mappings_to_load_on_current_rid, file);
        // std::cerr << "Load mappings\n";
        num_loaded_mappings_on_current_rid +=
            num_mappings_to_load_on_current_rid;
        num_mappings = num_mappings_to_load_on_current_rid;
      } else {
        // Move to next rid
        ++current_rid;
        if (current_rid < num_reference_sequences) {
          // std::cerr << "Load size\n";
          fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
          // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
          num_loaded_mappings_on_current_rid = 0;
        } else {
          break;
        }
      }
    }

    current_mapping_index = 0;
  }

  inline void Next(uint32_t num_reference_sequences) {
    ++current_mapping_index;
    if (current_mapping_index >= num_mappings) {
      LoadTempMappingBlock(num_reference_sequences);
    }
  }
};

template <>
inline void TempMappingFileHandle<PAFMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        break;
      }
    }
  }
  current_mapping_index = 0;
}

template <>
inline void TempMappingFileHandle<PairedPAFMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        break;
      }
    }
  }
  current_mapping_index = 0;
}

template <>
inline void TempMappingFileHandle<SAMMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        break;
      }
    }
  }
  current_mapping_index = 0;
}

template <>
inline void TempMappingFileHandle<PairsMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        break;
      }
    }
  }
  current_mapping_index = 0;
}

}  // namespace chromap

#endif  // TEMPMAPPING_H_
