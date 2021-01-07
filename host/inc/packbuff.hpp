#ifndef PACKBUFF_H_
#define PACKBUFF_H_

/*
  Class that handles a buffer to be packed with
  arrays of various sizes. The buffer shall be transfered
  to an accelerator device (eg. VE).
  Each array is allocated at an 8-byte aligned address.

  The procedure consists of the steps:
  - create a PackBuff() object
  - add arrays to it by calling .pack()
    - the data buffer grows as needed
    - store offset at the fixup address which will contain the VE address
  - get final size of buffer and allocate it on the VE
  - fixup pointer addresses by adding the buffer's VE address to each
    of the fixup addresses
  - transfer the data buffer to the VE

*/

#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdexcept>
#include <vector>

#define ALIGN8B(x) (((uint64_t)(x) + 7UL) & ~7UL)


class PackBuff {
private:
  char *data_;  // data of the packed buffer
  size_t size_; // current size of the buffer
  size_t size_alloc;  // allocated size of the buffer
  size_t chunksz;  // minimum data buffer growth granularity
  std::vector<uint64_t> fixup_; // addresses that need to be "fixed" by adding the VE
public:
  PackBuff(size_t n) : chunksz(n), size_(0) {
    data_ = (char *)malloc(chunksz);
    if (!data_)
      throw std::runtime_error("ERROR: failed to allocate initial PackBuff!");
    size_alloc = chunksz;
  }
  ~PackBuff() {
    if (data_)
      free(data_);
  }
  /*
   */
  void pack(void *buff, size_t len, uint64_t address) {
    size_t asize = ALIGN8B(size_);  // aligned pointer to buff destination
    // realloc if needed
    if (asize + len > size_alloc) {
      size_t new_size = ((asize + len + chunksz - 1) / chunksz ) * chunksz;
      char *new_data = (char *)realloc(data_, new_size);
      if (!new_data)
        throw std::runtime_error("ERROR: failed to realloc PackBuff!");

      // fix relocations which are inside the old buffer!
      for (auto it = fixup_.begin(); it != fixup_.end(); ++it) {
        auto addr = *it;
        if (addr >= (uint64_t)data_ && addr < (uint64_t)(data_ + size_)) {
          size_t offset = addr - (uint64_t)data_;
          *it = (uint64_t)(new_data + offset);
        }
      }
      // fix current argument's address if inside the relocated region
      if (address >= (uint64_t)data_ && address < (uint64_t)(data_ + size_)) {
        size_t offset = address - (uint64_t)data_;
        address = (uint64_t)(new_data + offset);
      }
      data_ = new_data;
      size_alloc = new_size;
    }
    if (buff) {
      // copy buffer to data buffer
      memcpy((void *)(data_ + asize), (const void *)buff, len);
    }
    // mark fixup address
    *((uint64_t *)address) = asize;
    fixup_.push_back(address);
    size_ = asize + len;
  }
  /*
   */
  size_t size() {
    return size_;
  }
  /*
   */
  void fixup(uint64_t VE_addr) {
    for (auto it = fixup_.begin(); it != fixup_.end(); ++it) {
      auto address = *it;
      *((uint64_t *)address) = *((uint64_t *)address) + VE_addr;
    }
  }
  /*
   */
  void *data() {
    return (void *)(data_);
  }
};

#endif // PACKBUFF_H_
