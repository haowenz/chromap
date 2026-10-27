#ifndef TEMPMAPPING_H_
#define TEMPMAPPING_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
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
  // The largest number of mappings on one ref seq in this file, recorded when
  // the file is written. A block never holds more than one ref seq's mappings,
  // so the loading buffer needs no more records than this.
  uint32_t max_num_mappings_on_a_rid = std::numeric_limits<uint32_t>::max();
  uint32_t current_rid;
  uint32_t current_mapping_index;
  uint32_t num_mappings_on_current_rid;
  uint32_t num_loaded_mappings_on_current_rid;
  // This vector only keep mappings on the same ref seq.
  std::vector<MappingRecord> mappings;
  // Record types of variable length (PairsMapping) are read through this
  // buffer, one fread per refill; read_buffer_position and read_buffer_end
  // delimit its unread bytes.
  std::vector<char> read_buffer;
  size_t read_buffer_position = 0;
  size_t read_buffer_end = 0;

  inline const MappingRecord& GetCurrentMapping() const {
    return mappings[current_mapping_index];
  }

  inline bool HasMappings() const { return num_mappings != 0; }

  inline void InitializeTempMappingLoading(uint32_t temp_mapping_block_size) {
    file = fopen(file_path.c_str(), "rb");
    if (file == NULL) {
      std::cerr << "Cannot open temporary file " << file_path << ". This may be caused by creating too many temporary files, please consider using command like \"ulimit -n 32768 -u 32768\" to increase the limit.\n" ;
    }
    assert(file != NULL);
    num_mappings = 0;
    block_size = temp_mapping_block_size;
    // Each block is loaded from one ref seq only, so capping the block size at
    // the largest ref seq's count loads the same records per block and avoids
    // allocating and zeroing a buffer far larger than the file.
    const uint32_t block_size_cap =
        max_num_mappings_on_a_rid > 0 ? max_num_mappings_on_a_rid : 1;
    if (block_size > block_size_cap) {
      block_size = block_size_cap;
    }
    current_rid = 0;
    current_mapping_index = 0;
    fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
    num_loaded_mappings_on_current_rid = 0;
    mappings.resize(block_size);
    read_buffer_position = 0;
    read_buffer_end = 0;
    // std::cerr << "Block size: " << block_size << ", initialize temp file " <<
    // file_path << "\n";
  }

  inline void FinalizeTempMappingLoading() { fclose(file); }

  // Returns a pointer to the next num_bytes unread bytes of the file without
  // consuming them, refilling read_buffer when fewer are buffered. Bytes past
  // the end of the file read as zero.
  inline const char* PeekBufferedBytes(size_t num_bytes) {
    if (read_buffer_end - read_buffer_position < num_bytes) {
      const size_t num_unread_bytes = read_buffer_end - read_buffer_position;
      if (num_unread_bytes > 0 && read_buffer_position > 0) {
        memmove(read_buffer.data(), read_buffer.data() + read_buffer_position,
                num_unread_bytes);
      }
      read_buffer_position = 0;
      read_buffer_end = num_unread_bytes;
      const size_t kReadBufferSize = (size_t)1 << 16;
      const size_t buffer_size =
          num_bytes > kReadBufferSize ? num_bytes : kReadBufferSize;
      if (read_buffer.size() < buffer_size) {
        read_buffer.resize(buffer_size);
      }
      read_buffer_end +=
          fread(read_buffer.data() + read_buffer_end, sizeof(char),
                read_buffer.size() - read_buffer_end, file);
      if (read_buffer_end < num_bytes) {
        memset(read_buffer.data() + read_buffer_end, 0,
               num_bytes - read_buffer_end);
        read_buffer_end = num_bytes;
      }
    }
    return read_buffer.data() + read_buffer_position;
  }

  // Reads the size_t count of mappings on the next ref seq from read_buffer.
  inline void LoadBufferedNumMappingsOnCurrentRid() {
    size_t num_mappings_on_rid = 0;
    memcpy(&num_mappings_on_rid, PeekBufferedBytes(sizeof(size_t)),
           sizeof(size_t));
    read_buffer_position += sizeof(size_t);
    num_mappings_on_current_rid = (uint32_t)num_mappings_on_rid;
  }

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
      // The bytes LoadFromFile reads, parsed from read_buffer.
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        const size_t num_record_bytes = PairsMapping::GetSerializedSize(
            PeekBufferedBytes(PairsMapping::kSerializedHeadSize));
        mappings[mi].DeserializeFrom(PeekBufferedBytes(num_record_bytes));
        read_buffer_position += num_record_bytes;
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
        LoadBufferedNumMappingsOnCurrentRid();
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
