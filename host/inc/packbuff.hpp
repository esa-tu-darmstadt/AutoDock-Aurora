/*

AutoDock-Aurora, a vectorized implementation of AutoDock 4.2 for the NEC SX-Aurora TSUBASA
Copyright (C) 2021 TU Darmstadt, Embedded Systems and Applications Group, Germany. All rights reserved.
For some of the code, Copyright (C) 2019 Computational Structural Biology Center, the Scripps Research Institute.

AutoDock is a Trade Mark of the Scripps Research Institute.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/

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

  Author: 2020-2021 Erich Focht, NEC
*/

#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>

#define ALIGN8B(x) (((uint64_t)(x) + 7UL) & ~7UL)


class PackBuff {
  struct reloc {
    uint64_t address;
    int64_t offset;
    reloc() : address(0), offset(0) {}
    reloc(uint64_t addr, uint64_t offs) : address(addr), offset(offs) {}
  };
private:
  char *data_;  // data of the packed buffer
  size_t size_; // current size of the buffer
  size_t size_alloc;  // allocated size of the buffer
  size_t chunksz;  // minimum data buffer growth granularity
  std::vector<reloc> fixup_; // addresses that need to be "fixed" by adding the VE
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
        if (it->address >= (uint64_t)data_ && it->address < (uint64_t)(data_ + size_)) {
          it->address = it->address - (uint64_t)data_ + (uint64_t)new_data;
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
    fixup_.push_back(reloc(address, asize));
    size_ = asize + len;
  }

  /*
    Return size of the data block.
   */
  size_t size() {
    return size_;
  }

  /*
    Save the relocation information and the data buffer into a file.
  */
  int save(const char* fname) {
    std::ofstream wf(fname, std::ios::out | std::ios::binary);
    if (!wf) {
      std::cerr << "Cannot open file!\n";
      return 1;
    }
    // count fixups with address inside the data block
    uint64_t data_begin = (uint64_t)this->data();
    uint64_t data_end = data_begin + this->size();
    int num_relocs = 0;
    struct reloc relocs[fixup_.size()];
    for (auto it = fixup_.begin(); it != fixup_.end(); ++it) {
      if (it->address >= data_begin && it->address < data_end) {
        relocs[num_relocs].address = it->address - data_begin;
        relocs[num_relocs].offset = it->offset;
        num_relocs++;
      }
    }
    // write length of data buffer
    wf.write((char *)&size_, sizeof(size_));
    // write number of relocations
    wf.write((char *)&num_relocs, sizeof(num_relocs));
    // write the relocations
    wf.write((char *)&relocs[0], sizeof(struct reloc) * num_relocs);
    // write data buffer
    wf.write(data_, this->size());
    wf.close();
    if(!wf.good()) {
      std::cerr << "Error occurred while saving!\n";
      return 1;
    }
  }

  /*
    Load the packed buffer and resolve relocations.
    This method returns a malloc'd data buffer that needs to be freed.
  */
  void *load(const char* fname) {
    std::ifstream rf(fname, std::ios::out | std::ios::binary);
    if(!rf) {
      std::cerr << "Cannot open 'load' file!\n";
      return nullptr;
    }

    size_t buff_size = 0;
    int num_relocs = 0;
    // read length of data buffer
    rf.read((char *)&buff_size, sizeof(buff_size));
    // allocate buffer
    void *buff = malloc(buff_size);
    if (buff == nullptr) {
      std::cerr << "Failed to malloc buffer of size " << buff_size
                << " bytes.\n";
      return nullptr;
    }
    
    // read number of relocations
    rf.read((char *)&num_relocs, sizeof(num_relocs));
    struct reloc relocs[num_relocs];

    // read the relocation offsets
    rf.read((char *)&relocs[0], sizeof(struct reloc) * num_relocs);

    // read the data buffer
    rf.read((char *)buff, buff_size);
    rf.close();
    if(!rf.good()) {
      std::cerr << "Error occurred while loading!\n";
      free(buff);
      return nullptr;
    }
    
    for (int i = 0; i < num_relocs; i++) {
      uint64_t address = (uint64_t)buff + relocs[i].address;
      *((uint64_t *)address) = relocs[i].offset + (uint64_t)buff;
    }
    return buff;
  }

  /*
    Fix relocations. The argument is the VE virtual address of the
    data buffer, i.e. the place where it shall be transfered. Relocations
    are added to this address.
   */
  void fixup(uint64_t VE_addr) {
    for (auto it = fixup_.begin(); it != fixup_.end(); ++it) {
      *((uint64_t *)it->address) = it->offset + VE_addr;
    }
  }
  /*
   */
  void *data() {
    return (void *)(data_);
  }
};

#endif // PACKBUFF_H_
