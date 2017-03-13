// (C) 1992-2016 Intel Corporation.                            
// Intel, the Intel logo, Intel, MegaCore, NIOS II, Quartus and TalkBack words    
// and logos are trademarks of Intel Corporation or its subsidiaries in the U.S.  
// and/or other countries. Other marks and brands may be claimed as the property  
// of others. See Trademarks on intel.com for full list of Intel trademarks or    
// the Trademarks & Brands Names Database (if Intel) or See www.Intel.com/legal (if Altera) 
// Your use of Intel Corporation's design tools, logic functions and other        
// software and tools, and its AMPP partner logic functions, and any output       
// files any of the foregoing (including device programming or simulation         
// files), and any associated documentation or information are expressly subject  
// to the terms and conditions of the Altera Program License Subscription         
// Agreement, Intel MegaCore Function License Agreement, or other applicable      
// license agreement, including, without limitation, that your use is for the     
// sole purpose of programming logic devices manufactured by Intel and sold by    
// Intel or its authorized distributors.  Please refer to the applicable          
// agreement for further details.                                                 
    


/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system
/////////////////////////////////////////////////////////////////
module Krnl_GA_system
(
   input logic clock,
   input logic clock2x,
   input logic resetn,
   // AVS avs_Krnl_Conform_cra
   input logic avs_Krnl_Conform_cra_enable,
   input logic avs_Krnl_Conform_cra_read,
   input logic avs_Krnl_Conform_cra_write,
   input logic [3:0] avs_Krnl_Conform_cra_address,
   input logic [63:0] avs_Krnl_Conform_cra_writedata,
   input logic [7:0] avs_Krnl_Conform_cra_byteenable,
   output logic [63:0] avs_Krnl_Conform_cra_readdata,
   output logic avs_Krnl_Conform_cra_readdatavalid,
   // AVS avs_Krnl_GA_cra
   input logic avs_Krnl_GA_cra_enable,
   input logic avs_Krnl_GA_cra_read,
   input logic avs_Krnl_GA_cra_write,
   input logic [4:0] avs_Krnl_GA_cra_address,
   input logic [63:0] avs_Krnl_GA_cra_writedata,
   input logic [7:0] avs_Krnl_GA_cra_byteenable,
   output logic [63:0] avs_Krnl_GA_cra_readdata,
   output logic avs_Krnl_GA_cra_readdatavalid,
   // AVS avs_Krnl_InterE_cra
   input logic avs_Krnl_InterE_cra_enable,
   input logic avs_Krnl_InterE_cra_read,
   input logic avs_Krnl_InterE_cra_write,
   input logic [4:0] avs_Krnl_InterE_cra_address,
   input logic [63:0] avs_Krnl_InterE_cra_writedata,
   input logic [7:0] avs_Krnl_InterE_cra_byteenable,
   output logic [63:0] avs_Krnl_InterE_cra_readdata,
   output logic avs_Krnl_InterE_cra_readdatavalid,
   // AVS avs_Krnl_IntraE_cra
   input logic avs_Krnl_IntraE_cra_enable,
   input logic avs_Krnl_IntraE_cra_read,
   input logic avs_Krnl_IntraE_cra_write,
   input logic [3:0] avs_Krnl_IntraE_cra_address,
   input logic [63:0] avs_Krnl_IntraE_cra_writedata,
   input logic [7:0] avs_Krnl_IntraE_cra_byteenable,
   output logic [63:0] avs_Krnl_IntraE_cra_readdata,
   output logic avs_Krnl_IntraE_cra_readdatavalid,
   // AVS avs_Krnl_Store_cra
   input logic avs_Krnl_Store_cra_enable,
   input logic avs_Krnl_Store_cra_read,
   input logic avs_Krnl_Store_cra_write,
   input logic [3:0] avs_Krnl_Store_cra_address,
   input logic [63:0] avs_Krnl_Store_cra_writedata,
   input logic [7:0] avs_Krnl_Store_cra_byteenable,
   output logic [63:0] avs_Krnl_Store_cra_readdata,
   output logic avs_Krnl_Store_cra_readdatavalid,
   output logic kernel_irq,
   // AVM avm_memgmem0_DDR_port_0_0_rw
   output logic avm_memgmem0_DDR_port_0_0_rw_enable,
   output logic avm_memgmem0_DDR_port_0_0_rw_read,
   output logic avm_memgmem0_DDR_port_0_0_rw_write,
   output logic [4:0] avm_memgmem0_DDR_port_0_0_rw_burstcount,
   output logic [30:0] avm_memgmem0_DDR_port_0_0_rw_address,
   output logic [511:0] avm_memgmem0_DDR_port_0_0_rw_writedata,
   output logic [63:0] avm_memgmem0_DDR_port_0_0_rw_byteenable,
   input logic avm_memgmem0_DDR_port_0_0_rw_waitrequest,
   input logic [511:0] avm_memgmem0_DDR_port_0_0_rw_readdata,
   input logic avm_memgmem0_DDR_port_0_0_rw_readdatavalid,
   input logic avm_memgmem0_DDR_port_0_0_rw_writeack
);
   genvar __i;
   logic Krnl_Conform_start;
   logic [0:0] Krnl_Conform_start_chain;
   logic [0:0] Krnl_Conform_start_kernel_copy;
   logic [0:0] Krnl_Conform_start_task_fd;
   logic [0:0] Krnl_Conform_start_finish_element;
   logic Krnl_Conform_finish;
   logic [0:0] Krnl_Conform_finish_kernel_copy;
   logic [0:0] Krnl_Conform_finish_chain;
   logic [31:0] Krnl_Conform_global_size [2:0];
   logic [31:0] Krnl_Conform_num_groups [2:0];
   logic [31:0] Krnl_Conform_local_size [2:0];
   logic [31:0] Krnl_Conform_global_offset [2:0];
   logic [31:0] Krnl_Conform_work_dim;
   logic [31:0] Krnl_Conform_wg_size;
   logic [0:0] Krnl_Conform_wg_disp_stall_in;
   logic [0:0] Krnl_Conform_wg_disp_valid_out;
   logic Krnl_Conform_wg_disp_start_out;
   logic [31:0] Krnl_Conform_wg_disp_group_id_out [2:0];
   logic [31:0] Krnl_Conform_wg_disp_global_id_base_out [2:0];
   logic Krnl_Conform_wg_disp_dispatched_all_groups;
   logic [31:0] Krnl_Conform_global_id [1][2:0];
   logic [31:0] Krnl_Conform_local_id [1][2:0];
   logic [31:0] Krnl_Conform_group_id [1][2:0];
   logic [0:0] Krnl_Conform_pending_write;
   logic [0:0] Krnl_Conform_lsu_active;
   logic [0:0] Krnl_Conform_valid_in;
   logic [0:0] Krnl_Conform_valid_out;
   logic [0:0] Krnl_Conform_stall_in;
   logic [0:0] Krnl_Conform_stall_out;
   logic Krnl_Conform_cra_pending_write;
   logic Krnl_Conform_cra_lsu_active;
   logic Krnl_Conform_cra_valid_in;
   logic [223:0] Krnl_Conform_kernel_arguments;
   logic Krnl_GA_start;
   logic [0:0] Krnl_GA_start_chain;
   logic [0:0] Krnl_GA_start_kernel_copy;
   logic [0:0] Krnl_GA_start_task_fd;
   logic [0:0] Krnl_GA_start_finish_element;
   logic Krnl_GA_finish;
   logic [0:0] Krnl_GA_finish_kernel_copy;
   logic [0:0] Krnl_GA_finish_chain;
   logic [31:0] Krnl_GA_global_size [2:0];
   logic [31:0] Krnl_GA_num_groups [2:0];
   logic [31:0] Krnl_GA_local_size [2:0];
   logic [31:0] Krnl_GA_global_offset [2:0];
   logic [31:0] Krnl_GA_work_dim;
   logic [31:0] Krnl_GA_wg_size;
   logic [0:0] Krnl_GA_wg_disp_stall_in;
   logic [0:0] Krnl_GA_wg_disp_valid_out;
   logic Krnl_GA_wg_disp_start_out;
   logic [31:0] Krnl_GA_wg_disp_group_id_out [2:0];
   logic [31:0] Krnl_GA_wg_disp_global_id_base_out [2:0];
   logic Krnl_GA_wg_disp_dispatched_all_groups;
   logic [31:0] Krnl_GA_global_id [1][2:0];
   logic [31:0] Krnl_GA_local_id [1][2:0];
   logic [31:0] Krnl_GA_group_id [1][2:0];
   logic [0:0] Krnl_GA_pending_write;
   logic [0:0] Krnl_GA_lsu_active;
   logic [0:0] Krnl_GA_valid_in;
   logic [0:0] Krnl_GA_valid_out;
   logic [0:0] Krnl_GA_stall_in;
   logic [0:0] Krnl_GA_stall_out;
   logic Krnl_GA_cra_pending_write;
   logic Krnl_GA_cra_lsu_active;
   logic Krnl_GA_cra_valid_in;
   logic [671:0] Krnl_GA_kernel_arguments;
   logic Krnl_InterE_start;
   logic [0:0] Krnl_InterE_start_chain;
   logic [0:0] Krnl_InterE_start_kernel_copy;
   logic [0:0] Krnl_InterE_start_task_fd;
   logic [0:0] Krnl_InterE_start_finish_element;
   logic Krnl_InterE_finish;
   logic [0:0] Krnl_InterE_finish_kernel_copy;
   logic [0:0] Krnl_InterE_finish_chain;
   logic [31:0] Krnl_InterE_global_size [2:0];
   logic [31:0] Krnl_InterE_num_groups [2:0];
   logic [31:0] Krnl_InterE_local_size [2:0];
   logic [31:0] Krnl_InterE_global_offset [2:0];
   logic [31:0] Krnl_InterE_work_dim;
   logic [31:0] Krnl_InterE_wg_size;
   logic [0:0] Krnl_InterE_wg_disp_stall_in;
   logic [0:0] Krnl_InterE_wg_disp_valid_out;
   logic Krnl_InterE_wg_disp_start_out;
   logic [31:0] Krnl_InterE_wg_disp_group_id_out [2:0];
   logic [31:0] Krnl_InterE_wg_disp_global_id_base_out [2:0];
   logic Krnl_InterE_wg_disp_dispatched_all_groups;
   logic [31:0] Krnl_InterE_global_id [1][2:0];
   logic [31:0] Krnl_InterE_local_id [1][2:0];
   logic [31:0] Krnl_InterE_group_id [1][2:0];
   logic [0:0] Krnl_InterE_pending_write;
   logic [0:0] Krnl_InterE_lsu_active;
   logic [0:0] Krnl_InterE_valid_in;
   logic [0:0] Krnl_InterE_valid_out;
   logic [0:0] Krnl_InterE_stall_in;
   logic [0:0] Krnl_InterE_stall_out;
   logic Krnl_InterE_cra_pending_write;
   logic Krnl_InterE_cra_lsu_active;
   logic Krnl_InterE_cra_valid_in;
   logic [287:0] Krnl_InterE_kernel_arguments;
   logic Krnl_IntraE_start;
   logic [0:0] Krnl_IntraE_start_chain;
   logic [0:0] Krnl_IntraE_start_kernel_copy;
   logic [0:0] Krnl_IntraE_start_task_fd;
   logic [0:0] Krnl_IntraE_start_finish_element;
   logic Krnl_IntraE_finish;
   logic [0:0] Krnl_IntraE_finish_kernel_copy;
   logic [0:0] Krnl_IntraE_finish_chain;
   logic [31:0] Krnl_IntraE_global_size [2:0];
   logic [31:0] Krnl_IntraE_num_groups [2:0];
   logic [31:0] Krnl_IntraE_local_size [2:0];
   logic [31:0] Krnl_IntraE_global_offset [2:0];
   logic [31:0] Krnl_IntraE_work_dim;
   logic [31:0] Krnl_IntraE_wg_size;
   logic [0:0] Krnl_IntraE_wg_disp_stall_in;
   logic [0:0] Krnl_IntraE_wg_disp_valid_out;
   logic Krnl_IntraE_wg_disp_start_out;
   logic [31:0] Krnl_IntraE_wg_disp_group_id_out [2:0];
   logic [31:0] Krnl_IntraE_wg_disp_global_id_base_out [2:0];
   logic Krnl_IntraE_wg_disp_dispatched_all_groups;
   logic [31:0] Krnl_IntraE_global_id [1][2:0];
   logic [31:0] Krnl_IntraE_local_id [1][2:0];
   logic [31:0] Krnl_IntraE_group_id [1][2:0];
   logic [0:0] Krnl_IntraE_pending_write;
   logic [0:0] Krnl_IntraE_lsu_active;
   logic [0:0] Krnl_IntraE_valid_in;
   logic [0:0] Krnl_IntraE_valid_out;
   logic [0:0] Krnl_IntraE_stall_in;
   logic [0:0] Krnl_IntraE_stall_out;
   logic Krnl_IntraE_cra_pending_write;
   logic Krnl_IntraE_cra_lsu_active;
   logic Krnl_IntraE_cra_valid_in;
   logic [223:0] Krnl_IntraE_kernel_arguments;
   logic Krnl_Store_start;
   logic [0:0] Krnl_Store_start_chain;
   logic [0:0] Krnl_Store_start_kernel_copy;
   logic [0:0] Krnl_Store_start_task_fd;
   logic [0:0] Krnl_Store_start_finish_element;
   logic Krnl_Store_finish;
   logic [0:0] Krnl_Store_finish_kernel_copy;
   logic [0:0] Krnl_Store_finish_chain;
   logic [31:0] Krnl_Store_global_size [2:0];
   logic [31:0] Krnl_Store_num_groups [2:0];
   logic [31:0] Krnl_Store_local_size [2:0];
   logic [31:0] Krnl_Store_global_offset [2:0];
   logic [31:0] Krnl_Store_work_dim;
   logic [31:0] Krnl_Store_wg_size;
   logic [0:0] Krnl_Store_wg_disp_stall_in;
   logic [0:0] Krnl_Store_wg_disp_valid_out;
   logic Krnl_Store_wg_disp_start_out;
   logic [31:0] Krnl_Store_wg_disp_group_id_out [2:0];
   logic [31:0] Krnl_Store_wg_disp_global_id_base_out [2:0];
   logic Krnl_Store_wg_disp_dispatched_all_groups;
   logic [31:0] Krnl_Store_global_id [1][2:0];
   logic [31:0] Krnl_Store_local_id [1][2:0];
   logic [31:0] Krnl_Store_group_id [1][2:0];
   logic [0:0] Krnl_Store_pending_write;
   logic [0:0] Krnl_Store_lsu_active;
   logic [0:0] Krnl_Store_valid_in;
   logic [0:0] Krnl_Store_valid_out;
   logic [0:0] Krnl_Store_stall_in;
   logic [0:0] Krnl_Store_stall_out;
   logic Krnl_Store_cra_pending_write;
   logic Krnl_Store_cra_lsu_active;
   logic Krnl_Store_cra_valid_in;
   logic [223:0] Krnl_Store_kernel_arguments;
   logic [4:0] kernel_irqs;
   logic avm_kernel_rd_enable [49];
   logic avm_kernel_rd_read [49];
   logic avm_kernel_rd_write [49];
   logic [4:0] avm_kernel_rd_burstcount [49];
   logic [30:0] avm_kernel_rd_address [49];
   logic [511:0] avm_kernel_rd_writedata [49];
   logic [63:0] avm_kernel_rd_byteenable [49];
   logic avm_kernel_rd_waitrequest [49];
   logic [511:0] avm_kernel_rd_readdata [49];
   logic avm_kernel_rd_readdatavalid [49];
   logic avm_kernel_rd_writeack [49];
   logic avm_kernel_wr_enable [10];
   logic avm_kernel_wr_read [10];
   logic avm_kernel_wr_write [10];
   logic [4:0] avm_kernel_wr_burstcount [10];
   logic [30:0] avm_kernel_wr_address [10];
   logic [511:0] avm_kernel_wr_writedata [10];
   logic [63:0] avm_kernel_wr_byteenable [10];
   logic avm_kernel_wr_waitrequest [10];
   logic [511:0] avm_kernel_wr_readdata [10];
   logic avm_kernel_wr_readdatavalid [10];
   logic avm_kernel_wr_writeack [10];
   logic ic_avm_enable [1];
   logic ic_avm_read [1];
   logic ic_avm_write [1];
   logic [4:0] ic_avm_burstcount [1];
   logic [30:0] ic_avm_address [1];
   logic [511:0] ic_avm_writedata [1];
   logic [63:0] ic_avm_byteenable [1];
   logic ic_avm_waitrequest [1];
   logic [511:0] ic_avm_readdata [1];
   logic ic_avm_readdatavalid [1];
   logic ic_avm_writeack [1];
   logic [31:0] avs_channel_id_chan_GA2Conf_genotype_fifosize;
   logic avm_channel_id_chan_GA2Conf_genotype_write_valid;
   logic avm_channel_id_chan_GA2Conf_genotype_write_ready;
   logic [31:0] avm_channel_id_chan_GA2Conf_genotype_write_data;
   logic avm_channel_id_chan_GA2Conf_genotype_read_valid;
   logic avm_channel_id_chan_GA2Conf_genotype_read_ready;
   logic [31:0] avm_channel_id_chan_GA2Conf_genotype_read_data;
   logic avs_channel_id_chan_GA2Conf_genotype_write_valid;
   logic avs_channel_id_chan_GA2Conf_genotype_write_ready;
   logic [31:0] avs_channel_id_chan_GA2Conf_genotype_write_data;
   logic avs_channel_id_chan_GA2Conf_genotype_read_valid;
   logic avs_channel_id_chan_GA2Conf_genotype_read_ready;
   logic [31:0] avs_channel_id_chan_GA2Conf_genotype_read_data;
   logic avs_channel_id_chan_GA2Conf_genotype_write_almostfull;
   logic avm_channel_id_chan_GA2Conf_genotype_write_almostfull;
   logic [31:0] avs_channel_id_chan_GA2Conf_active_fifosize;
   logic avm_channel_id_chan_GA2Conf_active_write_valid;
   logic avm_channel_id_chan_GA2Conf_active_write_ready;
   logic [7:0] avm_channel_id_chan_GA2Conf_active_write_data;
   logic avm_channel_id_chan_GA2Conf_active_read_valid;
   logic avm_channel_id_chan_GA2Conf_active_read_ready;
   logic [7:0] avm_channel_id_chan_GA2Conf_active_read_data;
   logic avs_channel_id_chan_GA2Conf_active_write_valid;
   logic avs_channel_id_chan_GA2Conf_active_write_ready;
   logic [7:0] avs_channel_id_chan_GA2Conf_active_write_data;
   logic avs_channel_id_chan_GA2Conf_active_read_valid;
   logic avs_channel_id_chan_GA2Conf_active_read_ready;
   logic [7:0] avs_channel_id_chan_GA2Conf_active_read_data;
   logic avs_channel_id_chan_GA2Conf_active_write_almostfull;
   logic avm_channel_id_chan_GA2Conf_active_write_almostfull;
   logic [31:0] avs_channel_id_chan_GA2Conf_cnt_fifosize;
   logic avm_channel_id_chan_GA2Conf_cnt_write_valid;
   logic avm_channel_id_chan_GA2Conf_cnt_write_ready;
   logic [31:0] avm_channel_id_chan_GA2Conf_cnt_write_data;
   logic avm_channel_id_chan_GA2Conf_cnt_read_valid;
   logic avm_channel_id_chan_GA2Conf_cnt_read_ready;
   logic [31:0] avm_channel_id_chan_GA2Conf_cnt_read_data;
   logic avs_channel_id_chan_GA2Conf_cnt_write_valid;
   logic avs_channel_id_chan_GA2Conf_cnt_write_ready;
   logic [31:0] avs_channel_id_chan_GA2Conf_cnt_write_data;
   logic avs_channel_id_chan_GA2Conf_cnt_read_valid;
   logic avs_channel_id_chan_GA2Conf_cnt_read_ready;
   logic [31:0] avs_channel_id_chan_GA2Conf_cnt_read_data;
   logic avs_channel_id_chan_GA2Conf_cnt_write_almostfull;
   logic avm_channel_id_chan_GA2Conf_cnt_write_almostfull;
   logic [31:0] avs_channel_id_chan_GA2Conf_mode_fifosize;
   logic avm_channel_id_chan_GA2Conf_mode_write_valid;
   logic avm_channel_id_chan_GA2Conf_mode_write_ready;
   logic [7:0] avm_channel_id_chan_GA2Conf_mode_write_data;
   logic avm_channel_id_chan_GA2Conf_mode_read_valid;
   logic avm_channel_id_chan_GA2Conf_mode_read_ready;
   logic [7:0] avm_channel_id_chan_GA2Conf_mode_read_data;
   logic avs_channel_id_chan_GA2Conf_mode_write_valid;
   logic avs_channel_id_chan_GA2Conf_mode_write_ready;
   logic [7:0] avs_channel_id_chan_GA2Conf_mode_write_data;
   logic avs_channel_id_chan_GA2Conf_mode_read_valid;
   logic avs_channel_id_chan_GA2Conf_mode_read_ready;
   logic [7:0] avs_channel_id_chan_GA2Conf_mode_read_data;
   logic avs_channel_id_chan_GA2Conf_mode_write_almostfull;
   logic avm_channel_id_chan_GA2Conf_mode_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intere_active_fifosize;
   logic avm_channel_id_chan_Conf2Intere_active_write_valid;
   logic avm_channel_id_chan_Conf2Intere_active_write_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intere_active_write_data;
   logic avm_channel_id_chan_Conf2Intere_active_read_valid;
   logic avm_channel_id_chan_Conf2Intere_active_read_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intere_active_read_data;
   logic avs_channel_id_chan_Conf2Intere_active_write_valid;
   logic avs_channel_id_chan_Conf2Intere_active_write_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intere_active_write_data;
   logic avs_channel_id_chan_Conf2Intere_active_read_valid;
   logic avs_channel_id_chan_Conf2Intere_active_read_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intere_active_read_data;
   logic avs_channel_id_chan_Conf2Intere_active_write_almostfull;
   logic avm_channel_id_chan_Conf2Intere_active_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intere_cnt_fifosize;
   logic avm_channel_id_chan_Conf2Intere_cnt_write_valid;
   logic avm_channel_id_chan_Conf2Intere_cnt_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_cnt_write_data;
   logic avm_channel_id_chan_Conf2Intere_cnt_read_valid;
   logic avm_channel_id_chan_Conf2Intere_cnt_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_cnt_read_data;
   logic avs_channel_id_chan_Conf2Intere_cnt_write_valid;
   logic avs_channel_id_chan_Conf2Intere_cnt_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_cnt_write_data;
   logic avs_channel_id_chan_Conf2Intere_cnt_read_valid;
   logic avs_channel_id_chan_Conf2Intere_cnt_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_cnt_read_data;
   logic avs_channel_id_chan_Conf2Intere_cnt_write_almostfull;
   logic avm_channel_id_chan_Conf2Intere_cnt_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intere_mode_fifosize;
   logic avm_channel_id_chan_Conf2Intere_mode_write_valid;
   logic avm_channel_id_chan_Conf2Intere_mode_write_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intere_mode_write_data;
   logic avm_channel_id_chan_Conf2Intere_mode_read_valid;
   logic avm_channel_id_chan_Conf2Intere_mode_read_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intere_mode_read_data;
   logic avs_channel_id_chan_Conf2Intere_mode_write_valid;
   logic avs_channel_id_chan_Conf2Intere_mode_write_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intere_mode_write_data;
   logic avs_channel_id_chan_Conf2Intere_mode_read_valid;
   logic avs_channel_id_chan_Conf2Intere_mode_read_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intere_mode_read_data;
   logic avs_channel_id_chan_Conf2Intere_mode_write_almostfull;
   logic avm_channel_id_chan_Conf2Intere_mode_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intere_x_fifosize;
   logic avm_channel_id_chan_Conf2Intere_x_write_valid;
   logic avm_channel_id_chan_Conf2Intere_x_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_x_write_data;
   logic avm_channel_id_chan_Conf2Intere_x_read_valid;
   logic avm_channel_id_chan_Conf2Intere_x_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_x_read_data;
   logic avs_channel_id_chan_Conf2Intere_x_write_valid;
   logic avs_channel_id_chan_Conf2Intere_x_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_x_write_data;
   logic avs_channel_id_chan_Conf2Intere_x_read_valid;
   logic avs_channel_id_chan_Conf2Intere_x_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_x_read_data;
   logic avs_channel_id_chan_Conf2Intere_x_write_almostfull;
   logic avm_channel_id_chan_Conf2Intere_x_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intere_y_fifosize;
   logic avm_channel_id_chan_Conf2Intere_y_write_valid;
   logic avm_channel_id_chan_Conf2Intere_y_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_y_write_data;
   logic avm_channel_id_chan_Conf2Intere_y_read_valid;
   logic avm_channel_id_chan_Conf2Intere_y_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_y_read_data;
   logic avs_channel_id_chan_Conf2Intere_y_write_valid;
   logic avs_channel_id_chan_Conf2Intere_y_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_y_write_data;
   logic avs_channel_id_chan_Conf2Intere_y_read_valid;
   logic avs_channel_id_chan_Conf2Intere_y_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_y_read_data;
   logic avs_channel_id_chan_Conf2Intere_y_write_almostfull;
   logic avm_channel_id_chan_Conf2Intere_y_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intere_z_fifosize;
   logic avm_channel_id_chan_Conf2Intere_z_write_valid;
   logic avm_channel_id_chan_Conf2Intere_z_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_z_write_data;
   logic avm_channel_id_chan_Conf2Intere_z_read_valid;
   logic avm_channel_id_chan_Conf2Intere_z_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intere_z_read_data;
   logic avs_channel_id_chan_Conf2Intere_z_write_valid;
   logic avs_channel_id_chan_Conf2Intere_z_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_z_write_data;
   logic avs_channel_id_chan_Conf2Intere_z_read_valid;
   logic avs_channel_id_chan_Conf2Intere_z_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intere_z_read_data;
   logic avs_channel_id_chan_Conf2Intere_z_write_almostfull;
   logic avm_channel_id_chan_Conf2Intere_z_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_active_fifosize;
   logic avm_channel_id_chan_Conf2Intrae_active_write_valid;
   logic avm_channel_id_chan_Conf2Intrae_active_write_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intrae_active_write_data;
   logic avm_channel_id_chan_Conf2Intrae_active_read_valid;
   logic avm_channel_id_chan_Conf2Intrae_active_read_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intrae_active_read_data;
   logic avs_channel_id_chan_Conf2Intrae_active_write_valid;
   logic avs_channel_id_chan_Conf2Intrae_active_write_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intrae_active_write_data;
   logic avs_channel_id_chan_Conf2Intrae_active_read_valid;
   logic avs_channel_id_chan_Conf2Intrae_active_read_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intrae_active_read_data;
   logic avs_channel_id_chan_Conf2Intrae_active_write_almostfull;
   logic avm_channel_id_chan_Conf2Intrae_active_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_cnt_fifosize;
   logic avm_channel_id_chan_Conf2Intrae_cnt_write_valid;
   logic avm_channel_id_chan_Conf2Intrae_cnt_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_cnt_write_data;
   logic avm_channel_id_chan_Conf2Intrae_cnt_read_valid;
   logic avm_channel_id_chan_Conf2Intrae_cnt_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_cnt_read_data;
   logic avs_channel_id_chan_Conf2Intrae_cnt_write_valid;
   logic avs_channel_id_chan_Conf2Intrae_cnt_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_cnt_write_data;
   logic avs_channel_id_chan_Conf2Intrae_cnt_read_valid;
   logic avs_channel_id_chan_Conf2Intrae_cnt_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_cnt_read_data;
   logic avs_channel_id_chan_Conf2Intrae_cnt_write_almostfull;
   logic avm_channel_id_chan_Conf2Intrae_cnt_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_mode_fifosize;
   logic avm_channel_id_chan_Conf2Intrae_mode_write_valid;
   logic avm_channel_id_chan_Conf2Intrae_mode_write_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intrae_mode_write_data;
   logic avm_channel_id_chan_Conf2Intrae_mode_read_valid;
   logic avm_channel_id_chan_Conf2Intrae_mode_read_ready;
   logic [7:0] avm_channel_id_chan_Conf2Intrae_mode_read_data;
   logic avs_channel_id_chan_Conf2Intrae_mode_write_valid;
   logic avs_channel_id_chan_Conf2Intrae_mode_write_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intrae_mode_write_data;
   logic avs_channel_id_chan_Conf2Intrae_mode_read_valid;
   logic avs_channel_id_chan_Conf2Intrae_mode_read_ready;
   logic [7:0] avs_channel_id_chan_Conf2Intrae_mode_read_data;
   logic avs_channel_id_chan_Conf2Intrae_mode_write_almostfull;
   logic avm_channel_id_chan_Conf2Intrae_mode_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_x_fifosize;
   logic avm_channel_id_chan_Conf2Intrae_x_write_valid;
   logic avm_channel_id_chan_Conf2Intrae_x_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_x_write_data;
   logic avm_channel_id_chan_Conf2Intrae_x_read_valid;
   logic avm_channel_id_chan_Conf2Intrae_x_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_x_read_data;
   logic avs_channel_id_chan_Conf2Intrae_x_write_valid;
   logic avs_channel_id_chan_Conf2Intrae_x_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_x_write_data;
   logic avs_channel_id_chan_Conf2Intrae_x_read_valid;
   logic avs_channel_id_chan_Conf2Intrae_x_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_x_read_data;
   logic avs_channel_id_chan_Conf2Intrae_x_write_almostfull;
   logic avm_channel_id_chan_Conf2Intrae_x_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_y_fifosize;
   logic avm_channel_id_chan_Conf2Intrae_y_write_valid;
   logic avm_channel_id_chan_Conf2Intrae_y_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_y_write_data;
   logic avm_channel_id_chan_Conf2Intrae_y_read_valid;
   logic avm_channel_id_chan_Conf2Intrae_y_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_y_read_data;
   logic avs_channel_id_chan_Conf2Intrae_y_write_valid;
   logic avs_channel_id_chan_Conf2Intrae_y_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_y_write_data;
   logic avs_channel_id_chan_Conf2Intrae_y_read_valid;
   logic avs_channel_id_chan_Conf2Intrae_y_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_y_read_data;
   logic avs_channel_id_chan_Conf2Intrae_y_write_almostfull;
   logic avm_channel_id_chan_Conf2Intrae_y_write_almostfull;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_z_fifosize;
   logic avm_channel_id_chan_Conf2Intrae_z_write_valid;
   logic avm_channel_id_chan_Conf2Intrae_z_write_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_z_write_data;
   logic avm_channel_id_chan_Conf2Intrae_z_read_valid;
   logic avm_channel_id_chan_Conf2Intrae_z_read_ready;
   logic [31:0] avm_channel_id_chan_Conf2Intrae_z_read_data;
   logic avs_channel_id_chan_Conf2Intrae_z_write_valid;
   logic avs_channel_id_chan_Conf2Intrae_z_write_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_z_write_data;
   logic avs_channel_id_chan_Conf2Intrae_z_read_valid;
   logic avs_channel_id_chan_Conf2Intrae_z_read_ready;
   logic [31:0] avs_channel_id_chan_Conf2Intrae_z_read_data;
   logic avs_channel_id_chan_Conf2Intrae_z_write_almostfull;
   logic avm_channel_id_chan_Conf2Intrae_z_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intere2Store_active_fifosize;
   logic avm_channel_id_chan_Intere2Store_active_write_valid;
   logic avm_channel_id_chan_Intere2Store_active_write_ready;
   logic [7:0] avm_channel_id_chan_Intere2Store_active_write_data;
   logic avm_channel_id_chan_Intere2Store_active_read_valid;
   logic avm_channel_id_chan_Intere2Store_active_read_ready;
   logic [7:0] avm_channel_id_chan_Intere2Store_active_read_data;
   logic avs_channel_id_chan_Intere2Store_active_write_valid;
   logic avs_channel_id_chan_Intere2Store_active_write_ready;
   logic [7:0] avs_channel_id_chan_Intere2Store_active_write_data;
   logic avs_channel_id_chan_Intere2Store_active_read_valid;
   logic avs_channel_id_chan_Intere2Store_active_read_ready;
   logic [7:0] avs_channel_id_chan_Intere2Store_active_read_data;
   logic avs_channel_id_chan_Intere2Store_active_write_almostfull;
   logic avm_channel_id_chan_Intere2Store_active_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intere2Store_cnt_fifosize;
   logic avm_channel_id_chan_Intere2Store_cnt_write_valid;
   logic avm_channel_id_chan_Intere2Store_cnt_write_ready;
   logic [31:0] avm_channel_id_chan_Intere2Store_cnt_write_data;
   logic avm_channel_id_chan_Intere2Store_cnt_read_valid;
   logic avm_channel_id_chan_Intere2Store_cnt_read_ready;
   logic [31:0] avm_channel_id_chan_Intere2Store_cnt_read_data;
   logic avs_channel_id_chan_Intere2Store_cnt_write_valid;
   logic avs_channel_id_chan_Intere2Store_cnt_write_ready;
   logic [31:0] avs_channel_id_chan_Intere2Store_cnt_write_data;
   logic avs_channel_id_chan_Intere2Store_cnt_read_valid;
   logic avs_channel_id_chan_Intere2Store_cnt_read_ready;
   logic [31:0] avs_channel_id_chan_Intere2Store_cnt_read_data;
   logic avs_channel_id_chan_Intere2Store_cnt_write_almostfull;
   logic avm_channel_id_chan_Intere2Store_cnt_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intere2Store_intere_fifosize;
   logic avm_channel_id_chan_Intere2Store_intere_write_valid;
   logic avm_channel_id_chan_Intere2Store_intere_write_ready;
   logic [31:0] avm_channel_id_chan_Intere2Store_intere_write_data;
   logic avm_channel_id_chan_Intere2Store_intere_read_valid;
   logic avm_channel_id_chan_Intere2Store_intere_read_ready;
   logic [31:0] avm_channel_id_chan_Intere2Store_intere_read_data;
   logic avs_channel_id_chan_Intere2Store_intere_write_valid;
   logic avs_channel_id_chan_Intere2Store_intere_write_ready;
   logic [31:0] avs_channel_id_chan_Intere2Store_intere_write_data;
   logic avs_channel_id_chan_Intere2Store_intere_read_valid;
   logic avs_channel_id_chan_Intere2Store_intere_read_ready;
   logic [31:0] avs_channel_id_chan_Intere2Store_intere_read_data;
   logic avs_channel_id_chan_Intere2Store_intere_write_almostfull;
   logic avm_channel_id_chan_Intere2Store_intere_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intere2Store_mode_fifosize;
   logic avm_channel_id_chan_Intere2Store_mode_write_valid;
   logic avm_channel_id_chan_Intere2Store_mode_write_ready;
   logic [7:0] avm_channel_id_chan_Intere2Store_mode_write_data;
   logic avm_channel_id_chan_Intere2Store_mode_read_valid;
   logic avm_channel_id_chan_Intere2Store_mode_read_ready;
   logic [7:0] avm_channel_id_chan_Intere2Store_mode_read_data;
   logic avs_channel_id_chan_Intere2Store_mode_write_valid;
   logic avs_channel_id_chan_Intere2Store_mode_write_ready;
   logic [7:0] avs_channel_id_chan_Intere2Store_mode_write_data;
   logic avs_channel_id_chan_Intere2Store_mode_read_valid;
   logic avs_channel_id_chan_Intere2Store_mode_read_ready;
   logic [7:0] avs_channel_id_chan_Intere2Store_mode_read_data;
   logic avs_channel_id_chan_Intere2Store_mode_write_almostfull;
   logic avm_channel_id_chan_Intere2Store_mode_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intrae2Store_active_fifosize;
   logic avm_channel_id_chan_Intrae2Store_active_write_valid;
   logic avm_channel_id_chan_Intrae2Store_active_write_ready;
   logic [7:0] avm_channel_id_chan_Intrae2Store_active_write_data;
   logic avm_channel_id_chan_Intrae2Store_active_read_valid;
   logic avm_channel_id_chan_Intrae2Store_active_read_ready;
   logic [7:0] avm_channel_id_chan_Intrae2Store_active_read_data;
   logic avs_channel_id_chan_Intrae2Store_active_write_valid;
   logic avs_channel_id_chan_Intrae2Store_active_write_ready;
   logic [7:0] avs_channel_id_chan_Intrae2Store_active_write_data;
   logic avs_channel_id_chan_Intrae2Store_active_read_valid;
   logic avs_channel_id_chan_Intrae2Store_active_read_ready;
   logic [7:0] avs_channel_id_chan_Intrae2Store_active_read_data;
   logic avs_channel_id_chan_Intrae2Store_active_write_almostfull;
   logic avm_channel_id_chan_Intrae2Store_active_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intrae2Store_cnt_fifosize;
   logic avm_channel_id_chan_Intrae2Store_cnt_write_valid;
   logic avm_channel_id_chan_Intrae2Store_cnt_write_ready;
   logic [31:0] avm_channel_id_chan_Intrae2Store_cnt_write_data;
   logic avm_channel_id_chan_Intrae2Store_cnt_read_valid;
   logic avm_channel_id_chan_Intrae2Store_cnt_read_ready;
   logic [31:0] avm_channel_id_chan_Intrae2Store_cnt_read_data;
   logic avs_channel_id_chan_Intrae2Store_cnt_write_valid;
   logic avs_channel_id_chan_Intrae2Store_cnt_write_ready;
   logic [31:0] avs_channel_id_chan_Intrae2Store_cnt_write_data;
   logic avs_channel_id_chan_Intrae2Store_cnt_read_valid;
   logic avs_channel_id_chan_Intrae2Store_cnt_read_ready;
   logic [31:0] avs_channel_id_chan_Intrae2Store_cnt_read_data;
   logic avs_channel_id_chan_Intrae2Store_cnt_write_almostfull;
   logic avm_channel_id_chan_Intrae2Store_cnt_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intrae2Store_intrae_fifosize;
   logic avm_channel_id_chan_Intrae2Store_intrae_write_valid;
   logic avm_channel_id_chan_Intrae2Store_intrae_write_ready;
   logic [31:0] avm_channel_id_chan_Intrae2Store_intrae_write_data;
   logic avm_channel_id_chan_Intrae2Store_intrae_read_valid;
   logic avm_channel_id_chan_Intrae2Store_intrae_read_ready;
   logic [31:0] avm_channel_id_chan_Intrae2Store_intrae_read_data;
   logic avs_channel_id_chan_Intrae2Store_intrae_write_valid;
   logic avs_channel_id_chan_Intrae2Store_intrae_write_ready;
   logic [31:0] avs_channel_id_chan_Intrae2Store_intrae_write_data;
   logic avs_channel_id_chan_Intrae2Store_intrae_read_valid;
   logic avs_channel_id_chan_Intrae2Store_intrae_read_ready;
   logic [31:0] avs_channel_id_chan_Intrae2Store_intrae_read_data;
   logic avs_channel_id_chan_Intrae2Store_intrae_write_almostfull;
   logic avm_channel_id_chan_Intrae2Store_intrae_write_almostfull;
   logic [31:0] avs_channel_id_chan_Intrae2Store_mode_fifosize;
   logic avm_channel_id_chan_Intrae2Store_mode_write_valid;
   logic avm_channel_id_chan_Intrae2Store_mode_write_ready;
   logic [7:0] avm_channel_id_chan_Intrae2Store_mode_write_data;
   logic avm_channel_id_chan_Intrae2Store_mode_read_valid;
   logic avm_channel_id_chan_Intrae2Store_mode_read_ready;
   logic [7:0] avm_channel_id_chan_Intrae2Store_mode_read_data;
   logic avs_channel_id_chan_Intrae2Store_mode_write_valid;
   logic avs_channel_id_chan_Intrae2Store_mode_write_ready;
   logic [7:0] avs_channel_id_chan_Intrae2Store_mode_write_data;
   logic avs_channel_id_chan_Intrae2Store_mode_read_valid;
   logic avs_channel_id_chan_Intrae2Store_mode_read_ready;
   logic [7:0] avs_channel_id_chan_Intrae2Store_mode_read_data;
   logic avs_channel_id_chan_Intrae2Store_mode_write_almostfull;
   logic avm_channel_id_chan_Intrae2Store_mode_write_almostfull;
   logic counter_reset_Krnl_Conform;
   logic [63:0] counter_init_Krnl_Conform;
   logic [31:0] counter_limit_Krnl_Conform;
   logic [31:0] counter_size_Krnl_Conform;
   logic counter_full_Krnl_Conform;
   logic counter_resetn_Krnl_Conform;
   logic [31:0] counter_increment_Krnl_Conform;
   logic counter_resultvalid_Krnl_Conform;
   logic counter_stall_Krnl_Conform;
   logic [63:0] counter_result_Krnl_Conform;
   logic counter_enable_Krnl_Conform;
   logic printf_avm_Krnl_Conform_enable [1];
   logic printf_avm_Krnl_Conform_read [1];
   logic printf_avm_Krnl_Conform_write [1];
   logic [5:0] printf_avm_Krnl_Conform_burstcount [1];
   logic [31:0] printf_avm_Krnl_Conform_address [1];
   logic [255:0] printf_avm_Krnl_Conform_writedata [1];
   logic [31:0] printf_avm_Krnl_Conform_byteenable [1];
   logic printf_avm_Krnl_Conform_waitrequest [1];
   logic [255:0] printf_avm_Krnl_Conform_readdata [1];
   logic printf_avm_Krnl_Conform_readdatavalid [1];
   logic counter_reset_Krnl_GA;
   logic [63:0] counter_init_Krnl_GA;
   logic [31:0] counter_limit_Krnl_GA;
   logic [31:0] counter_size_Krnl_GA;
   logic counter_full_Krnl_GA;
   logic counter_resetn_Krnl_GA;
   logic [31:0] counter_increment_Krnl_GA;
   logic counter_resultvalid_Krnl_GA;
   logic counter_stall_Krnl_GA;
   logic [63:0] counter_result_Krnl_GA;
   logic counter_enable_Krnl_GA;
   logic printf_avm_Krnl_GA_enable [1];
   logic printf_avm_Krnl_GA_read [1];
   logic printf_avm_Krnl_GA_write [1];
   logic [5:0] printf_avm_Krnl_GA_burstcount [1];
   logic [31:0] printf_avm_Krnl_GA_address [1];
   logic [255:0] printf_avm_Krnl_GA_writedata [1];
   logic [31:0] printf_avm_Krnl_GA_byteenable [1];
   logic printf_avm_Krnl_GA_waitrequest [1];
   logic [255:0] printf_avm_Krnl_GA_readdata [1];
   logic printf_avm_Krnl_GA_readdatavalid [1];
   logic counter_reset_Krnl_InterE;
   logic [63:0] counter_init_Krnl_InterE;
   logic [31:0] counter_limit_Krnl_InterE;
   logic [31:0] counter_size_Krnl_InterE;
   logic counter_full_Krnl_InterE;
   logic counter_resetn_Krnl_InterE;
   logic [31:0] counter_increment_Krnl_InterE;
   logic counter_resultvalid_Krnl_InterE;
   logic counter_stall_Krnl_InterE;
   logic [63:0] counter_result_Krnl_InterE;
   logic counter_enable_Krnl_InterE;
   logic printf_avm_Krnl_InterE_enable [1];
   logic printf_avm_Krnl_InterE_read [1];
   logic printf_avm_Krnl_InterE_write [1];
   logic [5:0] printf_avm_Krnl_InterE_burstcount [1];
   logic [31:0] printf_avm_Krnl_InterE_address [1];
   logic [255:0] printf_avm_Krnl_InterE_writedata [1];
   logic [31:0] printf_avm_Krnl_InterE_byteenable [1];
   logic printf_avm_Krnl_InterE_waitrequest [1];
   logic [255:0] printf_avm_Krnl_InterE_readdata [1];
   logic printf_avm_Krnl_InterE_readdatavalid [1];
   logic counter_reset_Krnl_IntraE;
   logic [63:0] counter_init_Krnl_IntraE;
   logic [31:0] counter_limit_Krnl_IntraE;
   logic [31:0] counter_size_Krnl_IntraE;
   logic counter_full_Krnl_IntraE;
   logic counter_resetn_Krnl_IntraE;
   logic [31:0] counter_increment_Krnl_IntraE;
   logic counter_resultvalid_Krnl_IntraE;
   logic counter_stall_Krnl_IntraE;
   logic [63:0] counter_result_Krnl_IntraE;
   logic counter_enable_Krnl_IntraE;
   logic printf_avm_Krnl_IntraE_enable [1];
   logic printf_avm_Krnl_IntraE_read [1];
   logic printf_avm_Krnl_IntraE_write [1];
   logic [5:0] printf_avm_Krnl_IntraE_burstcount [1];
   logic [31:0] printf_avm_Krnl_IntraE_address [1];
   logic [255:0] printf_avm_Krnl_IntraE_writedata [1];
   logic [31:0] printf_avm_Krnl_IntraE_byteenable [1];
   logic printf_avm_Krnl_IntraE_waitrequest [1];
   logic [255:0] printf_avm_Krnl_IntraE_readdata [1];
   logic printf_avm_Krnl_IntraE_readdatavalid [1];
   logic counter_reset_Krnl_Store;
   logic [63:0] counter_init_Krnl_Store;
   logic [31:0] counter_limit_Krnl_Store;
   logic [31:0] counter_size_Krnl_Store;
   logic counter_full_Krnl_Store;
   logic counter_resetn_Krnl_Store;
   logic [31:0] counter_increment_Krnl_Store;
   logic counter_resultvalid_Krnl_Store;
   logic counter_stall_Krnl_Store;
   logic [63:0] counter_result_Krnl_Store;
   logic counter_enable_Krnl_Store;
   logic printf_avm_Krnl_Store_enable [4];
   logic printf_avm_Krnl_Store_read [4];
   logic printf_avm_Krnl_Store_write [4];
   logic [5:0] printf_avm_Krnl_Store_burstcount [4];
   logic [31:0] printf_avm_Krnl_Store_address [4];
   logic [255:0] printf_avm_Krnl_Store_writedata [4];
   logic [31:0] printf_avm_Krnl_Store_byteenable [4];
   logic printf_avm_Krnl_Store_waitrequest [4];
   logic [255:0] printf_avm_Krnl_Store_readdata [4];
   logic printf_avm_Krnl_Store_readdatavalid [4];

   assign Krnl_Conform_start_chain[0] = Krnl_Conform_start;
   assign Krnl_Conform_finish_chain[0] = 1'b1;
   assign Krnl_Conform_cra_pending_write = |Krnl_Conform_pending_write;
   assign Krnl_Conform_cra_lsu_active = |Krnl_Conform_lsu_active;
   assign Krnl_Conform_cra_valid_in = |Krnl_Conform_valid_in;
   assign Krnl_Conform_stall_in = 0;
   // INST Krnl_Conform_workgroup_dispatcher of acl_work_group_dispatcher
   acl_work_group_dispatcher
   #(
      .WIDTH(32),
      .NUM_COPIES(1),
      .RUN_FOREVER(0)
   )
   Krnl_Conform_workgroup_dispatcher
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Conform_start),
      .num_groups(Krnl_Conform_num_groups),
      .local_size(Krnl_Conform_local_size),
      .stall_in(Krnl_Conform_wg_disp_stall_in),
      .valid_out(Krnl_Conform_wg_disp_valid_out),
      .group_id_out(Krnl_Conform_wg_disp_group_id_out),
      .global_id_base_out(Krnl_Conform_wg_disp_global_id_base_out),
      .start_out(Krnl_Conform_wg_disp_start_out),
      .dispatched_all_groups(Krnl_Conform_wg_disp_dispatched_all_groups)
   );

   // INST Krnl_Conform_finish_detector of acl_kernel_finish_detector
   acl_kernel_finish_detector
   #(
      .NUM_COPIES(1),
      .WG_SIZE_W(32),
      .GLOBAL_ID_W(32),
      .TESSELLATION_SIZE(19)
   )
   Krnl_Conform_finish_detector
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Conform_start),
      .wg_size(Krnl_Conform_wg_size),
      .wg_dispatch_valid_out(Krnl_Conform_wg_disp_valid_out),
      .wg_dispatch_stall_in(Krnl_Conform_wg_disp_stall_in),
      .dispatched_all_groups(Krnl_Conform_wg_disp_dispatched_all_groups),
      .kernel_copy_valid_out(Krnl_Conform_valid_out),
      .kernel_copy_stall_in(Krnl_Conform_stall_in),
      .pending_writes(Krnl_Conform_cra_pending_write),
      .finish(Krnl_Conform_finish)
   );

   // INST Krnl_Conform_cra_slave_inst of Krnl_Conform_function_cra_slave
   Krnl_Conform_function_cra_slave Krnl_Conform_cra_slave_inst
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Conform_start),
      .finish(Krnl_Conform_finish),
      .global_offset_0(Krnl_Conform_global_offset[0]),
      .global_offset_1(Krnl_Conform_global_offset[1]),
      .global_offset_2(Krnl_Conform_global_offset[2]),
      .work_dim(Krnl_Conform_work_dim),
      .has_a_lsu_active(Krnl_Conform_cra_lsu_active),
      .has_a_write_pending(Krnl_Conform_cra_pending_write),
      .valid_in(Krnl_Conform_cra_valid_in),
      .global_size_0(Krnl_Conform_global_size[0]),
      .global_size_1(Krnl_Conform_global_size[1]),
      .global_size_2(Krnl_Conform_global_size[2]),
      .num_groups_0(Krnl_Conform_num_groups[0]),
      .num_groups_1(Krnl_Conform_num_groups[1]),
      .num_groups_2(Krnl_Conform_num_groups[2]),
      .local_size_0(Krnl_Conform_local_size[0]),
      .local_size_1(Krnl_Conform_local_size[1]),
      .local_size_2(Krnl_Conform_local_size[2]),
      .workgroup_size(Krnl_Conform_wg_size),
      .kernel_arguments(Krnl_Conform_kernel_arguments),
      .cra_irq(kernel_irqs[0]),
      // AVS avs_cra
      .avs_cra_enable(avs_Krnl_Conform_cra_enable),
      .avs_cra_read(avs_Krnl_Conform_cra_read),
      .avs_cra_write(avs_Krnl_Conform_cra_write),
      .avs_cra_address(avs_Krnl_Conform_cra_address),
      .avs_cra_writedata(avs_Krnl_Conform_cra_writedata),
      .avs_cra_byteenable(avs_Krnl_Conform_cra_byteenable),
      .avs_cra_readdata(avs_Krnl_Conform_cra_readdata),
      .avs_cra_readdatavalid(avs_Krnl_Conform_cra_readdatavalid),
      .acl_counter_reset(counter_reset_Krnl_Conform),
      .acl_counter_init(counter_init_Krnl_Conform),
      .acl_counter_limit(counter_limit_Krnl_Conform),
      .acl_counter_size(counter_size_Krnl_Conform),
      .acl_counter_full(counter_full_Krnl_Conform)
   );

   // INST Krnl_Conform_id_iter_inst_0 of acl_id_iterator
   acl_id_iterator
   #(
      .WIDTH(32),
      .LOCAL_WIDTH_X(1),
      .LOCAL_WIDTH_Y(1),
      .LOCAL_WIDTH_Z(1),
      .ENABLE_TESSELLATION(1)
   )
   Krnl_Conform_id_iter_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Conform_wg_disp_start_out),
      .valid_in(Krnl_Conform_wg_disp_valid_out[0]),
      .stall_out(Krnl_Conform_wg_disp_stall_in[0]),
      .stall_in(Krnl_Conform_stall_out[0]),
      .valid_out(Krnl_Conform_valid_in[0]),
      .group_id_in(Krnl_Conform_wg_disp_group_id_out),
      .global_id_base_in(Krnl_Conform_wg_disp_global_id_base_out),
      .local_size(Krnl_Conform_local_size),
      .global_size(Krnl_Conform_global_size),
      .local_id(Krnl_Conform_local_id[0]),
      .global_id(Krnl_Conform_global_id[0]),
      .group_id(Krnl_Conform_group_id[0])
   );

   // INST Krnl_Conform_inst_0 of Krnl_Conform_top_wrapper_0
   Krnl_Conform_top_wrapper_0 Krnl_Conform_inst_0
   (
      .start(Krnl_Conform_start_kernel_copy[0]),
      .kernel_arguments(Krnl_Conform_kernel_arguments),
      .work_dim(Krnl_Conform_work_dim),
      .global_offset(Krnl_Conform_global_offset),
      .kernel_valid_out(Krnl_Conform_valid_out[0]),
      .has_a_write_pending(Krnl_Conform_pending_write[0]),
      .has_a_lsu_active(Krnl_Conform_lsu_active[0]),
      .global_id(Krnl_Conform_global_id[0]),
      .local_id(Krnl_Conform_local_id[0]),
      .group_id(Krnl_Conform_group_id[0]),
      .global_size(Krnl_Conform_global_size),
      .local_size(Krnl_Conform_local_size),
      .num_groups(Krnl_Conform_num_groups),
      .workgroup_size(Krnl_Conform_wg_size),
      .kernel_stall_out(Krnl_Conform_stall_out[0]),
      .kernel_valid_in(Krnl_Conform_valid_in[0]),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable(avm_kernel_rd_enable[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read(avm_kernel_rd_read[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write(avm_kernel_rd_write[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount(avm_kernel_rd_burstcount[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address(avm_kernel_rd_address[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata(avm_kernel_rd_writedata[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable(avm_kernel_rd_byteenable[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest(avm_kernel_rd_waitrequest[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata(avm_kernel_rd_readdata[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid(avm_kernel_rd_readdatavalid[0]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack(avm_kernel_rd_writeack[0]),
      // AVM avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_enable(avm_kernel_rd_enable[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_read(avm_kernel_rd_read[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_write(avm_kernel_rd_write[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_burstcount(avm_kernel_rd_burstcount[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_address(avm_kernel_rd_address[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writedata(avm_kernel_rd_writedata[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_byteenable(avm_kernel_rd_byteenable[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_waitrequest(avm_kernel_rd_waitrequest[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdata(avm_kernel_rd_readdata[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid(avm_kernel_rd_readdatavalid[1]),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writeack(avm_kernel_rd_writeack[1]),
      // AVM avm_local_bb4_st__inst0
      .avm_local_bb4_st__inst0_enable(avm_kernel_wr_enable[0]),
      .avm_local_bb4_st__inst0_read(avm_kernel_wr_read[0]),
      .avm_local_bb4_st__inst0_write(avm_kernel_wr_write[0]),
      .avm_local_bb4_st__inst0_burstcount(avm_kernel_wr_burstcount[0]),
      .avm_local_bb4_st__inst0_address(avm_kernel_wr_address[0]),
      .avm_local_bb4_st__inst0_writedata(avm_kernel_wr_writedata[0]),
      .avm_local_bb4_st__inst0_byteenable(avm_kernel_wr_byteenable[0]),
      .avm_local_bb4_st__inst0_waitrequest(avm_kernel_wr_waitrequest[0]),
      .avm_local_bb4_st__inst0_readdata(avm_kernel_wr_readdata[0]),
      .avm_local_bb4_st__inst0_readdatavalid(avm_kernel_wr_readdatavalid[0]),
      .avm_local_bb4_st__inst0_writeack(avm_kernel_wr_writeack[0]),
      // AVM avm_local_bb5_ld__inst0
      .avm_local_bb5_ld__inst0_enable(avm_kernel_rd_enable[2]),
      .avm_local_bb5_ld__inst0_read(avm_kernel_rd_read[2]),
      .avm_local_bb5_ld__inst0_write(avm_kernel_rd_write[2]),
      .avm_local_bb5_ld__inst0_burstcount(avm_kernel_rd_burstcount[2]),
      .avm_local_bb5_ld__inst0_address(avm_kernel_rd_address[2]),
      .avm_local_bb5_ld__inst0_writedata(avm_kernel_rd_writedata[2]),
      .avm_local_bb5_ld__inst0_byteenable(avm_kernel_rd_byteenable[2]),
      .avm_local_bb5_ld__inst0_waitrequest(avm_kernel_rd_waitrequest[2]),
      .avm_local_bb5_ld__inst0_readdata(avm_kernel_rd_readdata[2]),
      .avm_local_bb5_ld__inst0_readdatavalid(avm_kernel_rd_readdatavalid[2]),
      .avm_local_bb5_ld__inst0_writeack(avm_kernel_rd_writeack[2]),
      // AVM avm_local_bb5_ld__u15_inst0
      .avm_local_bb5_ld__u15_inst0_enable(avm_kernel_rd_enable[3]),
      .avm_local_bb5_ld__u15_inst0_read(avm_kernel_rd_read[3]),
      .avm_local_bb5_ld__u15_inst0_write(avm_kernel_rd_write[3]),
      .avm_local_bb5_ld__u15_inst0_burstcount(avm_kernel_rd_burstcount[3]),
      .avm_local_bb5_ld__u15_inst0_address(avm_kernel_rd_address[3]),
      .avm_local_bb5_ld__u15_inst0_writedata(avm_kernel_rd_writedata[3]),
      .avm_local_bb5_ld__u15_inst0_byteenable(avm_kernel_rd_byteenable[3]),
      .avm_local_bb5_ld__u15_inst0_waitrequest(avm_kernel_rd_waitrequest[3]),
      .avm_local_bb5_ld__u15_inst0_readdata(avm_kernel_rd_readdata[3]),
      .avm_local_bb5_ld__u15_inst0_readdatavalid(avm_kernel_rd_readdatavalid[3]),
      .avm_local_bb5_ld__u15_inst0_writeack(avm_kernel_rd_writeack[3]),
      // AVM avm_local_bb5_ld__u16_inst0
      .avm_local_bb5_ld__u16_inst0_enable(avm_kernel_rd_enable[4]),
      .avm_local_bb5_ld__u16_inst0_read(avm_kernel_rd_read[4]),
      .avm_local_bb5_ld__u16_inst0_write(avm_kernel_rd_write[4]),
      .avm_local_bb5_ld__u16_inst0_burstcount(avm_kernel_rd_burstcount[4]),
      .avm_local_bb5_ld__u16_inst0_address(avm_kernel_rd_address[4]),
      .avm_local_bb5_ld__u16_inst0_writedata(avm_kernel_rd_writedata[4]),
      .avm_local_bb5_ld__u16_inst0_byteenable(avm_kernel_rd_byteenable[4]),
      .avm_local_bb5_ld__u16_inst0_waitrequest(avm_kernel_rd_waitrequest[4]),
      .avm_local_bb5_ld__u16_inst0_readdata(avm_kernel_rd_readdata[4]),
      .avm_local_bb5_ld__u16_inst0_readdatavalid(avm_kernel_rd_readdatavalid[4]),
      .avm_local_bb5_ld__u16_inst0_writeack(avm_kernel_rd_writeack[4]),
      // AVM avm_local_bb5_ld__u17_inst0
      .avm_local_bb5_ld__u17_inst0_enable(avm_kernel_rd_enable[5]),
      .avm_local_bb5_ld__u17_inst0_read(avm_kernel_rd_read[5]),
      .avm_local_bb5_ld__u17_inst0_write(avm_kernel_rd_write[5]),
      .avm_local_bb5_ld__u17_inst0_burstcount(avm_kernel_rd_burstcount[5]),
      .avm_local_bb5_ld__u17_inst0_address(avm_kernel_rd_address[5]),
      .avm_local_bb5_ld__u17_inst0_writedata(avm_kernel_rd_writedata[5]),
      .avm_local_bb5_ld__u17_inst0_byteenable(avm_kernel_rd_byteenable[5]),
      .avm_local_bb5_ld__u17_inst0_waitrequest(avm_kernel_rd_waitrequest[5]),
      .avm_local_bb5_ld__u17_inst0_readdata(avm_kernel_rd_readdata[5]),
      .avm_local_bb5_ld__u17_inst0_readdatavalid(avm_kernel_rd_readdatavalid[5]),
      .avm_local_bb5_ld__u17_inst0_writeack(avm_kernel_rd_writeack[5]),
      // AVM avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_enable(avm_kernel_rd_enable[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_read(avm_kernel_rd_read[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_write(avm_kernel_rd_write[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_burstcount(avm_kernel_rd_burstcount[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_address(avm_kernel_rd_address[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writedata(avm_kernel_rd_writedata[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_byteenable(avm_kernel_rd_byteenable[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_waitrequest(avm_kernel_rd_waitrequest[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdata(avm_kernel_rd_readdata[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdatavalid(avm_kernel_rd_readdatavalid[6]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writeack(avm_kernel_rd_writeack[6]),
      // AVM avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_enable(avm_kernel_rd_enable[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_read(avm_kernel_rd_read[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_write(avm_kernel_rd_write[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_burstcount(avm_kernel_rd_burstcount[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_address(avm_kernel_rd_address[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writedata(avm_kernel_rd_writedata[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_byteenable(avm_kernel_rd_byteenable[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_waitrequest(avm_kernel_rd_waitrequest[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdata(avm_kernel_rd_readdata[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdatavalid(avm_kernel_rd_readdatavalid[7]),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writeack(avm_kernel_rd_writeack[7]),
      // AVST avm_channel_id_chan_GA2Conf_genotype_read
      .avm_channel_id_chan_GA2Conf_genotype_read_valid(avm_channel_id_chan_GA2Conf_genotype_read_valid),
      .avm_channel_id_chan_GA2Conf_genotype_read_ready(avm_channel_id_chan_GA2Conf_genotype_read_ready),
      .avm_channel_id_chan_GA2Conf_genotype_read_data(avm_channel_id_chan_GA2Conf_genotype_read_data),
      // AVST avm_channel_id_chan_GA2Conf_active_read
      .avm_channel_id_chan_GA2Conf_active_read_valid(avm_channel_id_chan_GA2Conf_active_read_valid),
      .avm_channel_id_chan_GA2Conf_active_read_ready(avm_channel_id_chan_GA2Conf_active_read_ready),
      .avm_channel_id_chan_GA2Conf_active_read_data(avm_channel_id_chan_GA2Conf_active_read_data),
      // AVST avm_channel_id_chan_GA2Conf_cnt_read
      .avm_channel_id_chan_GA2Conf_cnt_read_valid(avm_channel_id_chan_GA2Conf_cnt_read_valid),
      .avm_channel_id_chan_GA2Conf_cnt_read_ready(avm_channel_id_chan_GA2Conf_cnt_read_ready),
      .avm_channel_id_chan_GA2Conf_cnt_read_data(avm_channel_id_chan_GA2Conf_cnt_read_data),
      // AVST avm_channel_id_chan_GA2Conf_mode_read
      .avm_channel_id_chan_GA2Conf_mode_read_valid(avm_channel_id_chan_GA2Conf_mode_read_valid),
      .avm_channel_id_chan_GA2Conf_mode_read_ready(avm_channel_id_chan_GA2Conf_mode_read_ready),
      .avm_channel_id_chan_GA2Conf_mode_read_data(avm_channel_id_chan_GA2Conf_mode_read_data),
      // AVST avm_channel_id_chan_Conf2Intere_active_write
      .avm_channel_id_chan_Conf2Intere_active_write_valid(avm_channel_id_chan_Conf2Intere_active_write_valid),
      .avm_channel_id_chan_Conf2Intere_active_write_ready(avm_channel_id_chan_Conf2Intere_active_write_ready),
      .avm_channel_id_chan_Conf2Intere_active_write_data(avm_channel_id_chan_Conf2Intere_active_write_data),
      .avm_channel_id_chan_Conf2Intere_active_write_almostfull(avm_channel_id_chan_Conf2Intere_active_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intere_cnt_write
      .avm_channel_id_chan_Conf2Intere_cnt_write_valid(avm_channel_id_chan_Conf2Intere_cnt_write_valid),
      .avm_channel_id_chan_Conf2Intere_cnt_write_ready(avm_channel_id_chan_Conf2Intere_cnt_write_ready),
      .avm_channel_id_chan_Conf2Intere_cnt_write_data(avm_channel_id_chan_Conf2Intere_cnt_write_data),
      .avm_channel_id_chan_Conf2Intere_cnt_write_almostfull(avm_channel_id_chan_Conf2Intere_cnt_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intere_mode_write
      .avm_channel_id_chan_Conf2Intere_mode_write_valid(avm_channel_id_chan_Conf2Intere_mode_write_valid),
      .avm_channel_id_chan_Conf2Intere_mode_write_ready(avm_channel_id_chan_Conf2Intere_mode_write_ready),
      .avm_channel_id_chan_Conf2Intere_mode_write_data(avm_channel_id_chan_Conf2Intere_mode_write_data),
      .avm_channel_id_chan_Conf2Intere_mode_write_almostfull(avm_channel_id_chan_Conf2Intere_mode_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intere_x_write
      .avm_channel_id_chan_Conf2Intere_x_write_valid(avm_channel_id_chan_Conf2Intere_x_write_valid),
      .avm_channel_id_chan_Conf2Intere_x_write_ready(avm_channel_id_chan_Conf2Intere_x_write_ready),
      .avm_channel_id_chan_Conf2Intere_x_write_data(avm_channel_id_chan_Conf2Intere_x_write_data),
      .avm_channel_id_chan_Conf2Intere_x_write_almostfull(avm_channel_id_chan_Conf2Intere_x_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intere_y_write
      .avm_channel_id_chan_Conf2Intere_y_write_valid(avm_channel_id_chan_Conf2Intere_y_write_valid),
      .avm_channel_id_chan_Conf2Intere_y_write_ready(avm_channel_id_chan_Conf2Intere_y_write_ready),
      .avm_channel_id_chan_Conf2Intere_y_write_data(avm_channel_id_chan_Conf2Intere_y_write_data),
      .avm_channel_id_chan_Conf2Intere_y_write_almostfull(avm_channel_id_chan_Conf2Intere_y_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intere_z_write
      .avm_channel_id_chan_Conf2Intere_z_write_valid(avm_channel_id_chan_Conf2Intere_z_write_valid),
      .avm_channel_id_chan_Conf2Intere_z_write_ready(avm_channel_id_chan_Conf2Intere_z_write_ready),
      .avm_channel_id_chan_Conf2Intere_z_write_data(avm_channel_id_chan_Conf2Intere_z_write_data),
      .avm_channel_id_chan_Conf2Intere_z_write_almostfull(avm_channel_id_chan_Conf2Intere_z_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intrae_active_write
      .avm_channel_id_chan_Conf2Intrae_active_write_valid(avm_channel_id_chan_Conf2Intrae_active_write_valid),
      .avm_channel_id_chan_Conf2Intrae_active_write_ready(avm_channel_id_chan_Conf2Intrae_active_write_ready),
      .avm_channel_id_chan_Conf2Intrae_active_write_data(avm_channel_id_chan_Conf2Intrae_active_write_data),
      .avm_channel_id_chan_Conf2Intrae_active_write_almostfull(avm_channel_id_chan_Conf2Intrae_active_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intrae_cnt_write
      .avm_channel_id_chan_Conf2Intrae_cnt_write_valid(avm_channel_id_chan_Conf2Intrae_cnt_write_valid),
      .avm_channel_id_chan_Conf2Intrae_cnt_write_ready(avm_channel_id_chan_Conf2Intrae_cnt_write_ready),
      .avm_channel_id_chan_Conf2Intrae_cnt_write_data(avm_channel_id_chan_Conf2Intrae_cnt_write_data),
      .avm_channel_id_chan_Conf2Intrae_cnt_write_almostfull(avm_channel_id_chan_Conf2Intrae_cnt_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intrae_mode_write
      .avm_channel_id_chan_Conf2Intrae_mode_write_valid(avm_channel_id_chan_Conf2Intrae_mode_write_valid),
      .avm_channel_id_chan_Conf2Intrae_mode_write_ready(avm_channel_id_chan_Conf2Intrae_mode_write_ready),
      .avm_channel_id_chan_Conf2Intrae_mode_write_data(avm_channel_id_chan_Conf2Intrae_mode_write_data),
      .avm_channel_id_chan_Conf2Intrae_mode_write_almostfull(avm_channel_id_chan_Conf2Intrae_mode_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intrae_x_write
      .avm_channel_id_chan_Conf2Intrae_x_write_valid(avm_channel_id_chan_Conf2Intrae_x_write_valid),
      .avm_channel_id_chan_Conf2Intrae_x_write_ready(avm_channel_id_chan_Conf2Intrae_x_write_ready),
      .avm_channel_id_chan_Conf2Intrae_x_write_data(avm_channel_id_chan_Conf2Intrae_x_write_data),
      .avm_channel_id_chan_Conf2Intrae_x_write_almostfull(avm_channel_id_chan_Conf2Intrae_x_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intrae_y_write
      .avm_channel_id_chan_Conf2Intrae_y_write_valid(avm_channel_id_chan_Conf2Intrae_y_write_valid),
      .avm_channel_id_chan_Conf2Intrae_y_write_ready(avm_channel_id_chan_Conf2Intrae_y_write_ready),
      .avm_channel_id_chan_Conf2Intrae_y_write_data(avm_channel_id_chan_Conf2Intrae_y_write_data),
      .avm_channel_id_chan_Conf2Intrae_y_write_almostfull(avm_channel_id_chan_Conf2Intrae_y_write_almostfull),
      // AVST avm_channel_id_chan_Conf2Intrae_z_write
      .avm_channel_id_chan_Conf2Intrae_z_write_valid(avm_channel_id_chan_Conf2Intrae_z_write_valid),
      .avm_channel_id_chan_Conf2Intrae_z_write_ready(avm_channel_id_chan_Conf2Intrae_z_write_ready),
      .avm_channel_id_chan_Conf2Intrae_z_write_data(avm_channel_id_chan_Conf2Intrae_z_write_data),
      .avm_channel_id_chan_Conf2Intrae_z_write_almostfull(avm_channel_id_chan_Conf2Intrae_z_write_almostfull),
      // AVM p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_Conform_enable[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_Conform_read[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_Conform_write[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_Conform_burstcount[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_Conform_address[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_Conform_writedata[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_Conform_byteenable[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_Conform_waitrequest[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_Conform_readdata[0]),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_Conform_readdatavalid[0])
   );

   // INST Krnl_Conform_start_elem_inst_0 of acl_start_signal_chain_element
   acl_start_signal_chain_element Krnl_Conform_start_elem_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start_in(Krnl_Conform_start_chain[0]),
      .start_kernel(Krnl_Conform_start_kernel_copy[0]),
      .start_finish_detector(Krnl_Conform_start_task_fd[0]),
      .start_finish_chain_element(Krnl_Conform_start_finish_element[0]),
      .start_chain()
   );

   assign Krnl_GA_start_chain[0] = Krnl_GA_start;
   assign Krnl_GA_finish_chain[0] = 1'b1;
   assign Krnl_GA_cra_pending_write = |Krnl_GA_pending_write;
   assign Krnl_GA_cra_lsu_active = |Krnl_GA_lsu_active;
   assign Krnl_GA_cra_valid_in = |Krnl_GA_valid_in;
   assign Krnl_GA_stall_in = 0;
   // INST Krnl_GA_workgroup_dispatcher of acl_work_group_dispatcher
   acl_work_group_dispatcher
   #(
      .WIDTH(32),
      .NUM_COPIES(1),
      .RUN_FOREVER(0)
   )
   Krnl_GA_workgroup_dispatcher
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_GA_start),
      .num_groups(Krnl_GA_num_groups),
      .local_size(Krnl_GA_local_size),
      .stall_in(Krnl_GA_wg_disp_stall_in),
      .valid_out(Krnl_GA_wg_disp_valid_out),
      .group_id_out(Krnl_GA_wg_disp_group_id_out),
      .global_id_base_out(Krnl_GA_wg_disp_global_id_base_out),
      .start_out(Krnl_GA_wg_disp_start_out),
      .dispatched_all_groups(Krnl_GA_wg_disp_dispatched_all_groups)
   );

   // INST Krnl_GA_finish_detector of acl_kernel_finish_detector
   acl_kernel_finish_detector
   #(
      .NUM_COPIES(1),
      .WG_SIZE_W(32),
      .GLOBAL_ID_W(32),
      .TESSELLATION_SIZE(19)
   )
   Krnl_GA_finish_detector
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_GA_start),
      .wg_size(Krnl_GA_wg_size),
      .wg_dispatch_valid_out(Krnl_GA_wg_disp_valid_out),
      .wg_dispatch_stall_in(Krnl_GA_wg_disp_stall_in),
      .dispatched_all_groups(Krnl_GA_wg_disp_dispatched_all_groups),
      .kernel_copy_valid_out(Krnl_GA_valid_out),
      .kernel_copy_stall_in(Krnl_GA_stall_in),
      .pending_writes(Krnl_GA_cra_pending_write),
      .finish(Krnl_GA_finish)
   );

   // INST Krnl_GA_cra_slave_inst of Krnl_GA_function_cra_slave
   Krnl_GA_function_cra_slave Krnl_GA_cra_slave_inst
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_GA_start),
      .finish(Krnl_GA_finish),
      .global_offset_0(Krnl_GA_global_offset[0]),
      .global_offset_1(Krnl_GA_global_offset[1]),
      .global_offset_2(Krnl_GA_global_offset[2]),
      .work_dim(Krnl_GA_work_dim),
      .has_a_lsu_active(Krnl_GA_cra_lsu_active),
      .has_a_write_pending(Krnl_GA_cra_pending_write),
      .valid_in(Krnl_GA_cra_valid_in),
      .global_size_0(Krnl_GA_global_size[0]),
      .global_size_1(Krnl_GA_global_size[1]),
      .global_size_2(Krnl_GA_global_size[2]),
      .num_groups_0(Krnl_GA_num_groups[0]),
      .num_groups_1(Krnl_GA_num_groups[1]),
      .num_groups_2(Krnl_GA_num_groups[2]),
      .local_size_0(Krnl_GA_local_size[0]),
      .local_size_1(Krnl_GA_local_size[1]),
      .local_size_2(Krnl_GA_local_size[2]),
      .workgroup_size(Krnl_GA_wg_size),
      .kernel_arguments(Krnl_GA_kernel_arguments),
      .cra_irq(kernel_irqs[1]),
      // AVS avs_cra
      .avs_cra_enable(avs_Krnl_GA_cra_enable),
      .avs_cra_read(avs_Krnl_GA_cra_read),
      .avs_cra_write(avs_Krnl_GA_cra_write),
      .avs_cra_address(avs_Krnl_GA_cra_address),
      .avs_cra_writedata(avs_Krnl_GA_cra_writedata),
      .avs_cra_byteenable(avs_Krnl_GA_cra_byteenable),
      .avs_cra_readdata(avs_Krnl_GA_cra_readdata),
      .avs_cra_readdatavalid(avs_Krnl_GA_cra_readdatavalid),
      .acl_counter_reset(counter_reset_Krnl_GA),
      .acl_counter_init(counter_init_Krnl_GA),
      .acl_counter_limit(counter_limit_Krnl_GA),
      .acl_counter_size(counter_size_Krnl_GA),
      .acl_counter_full(counter_full_Krnl_GA)
   );

   // INST Krnl_GA_id_iter_inst_0 of acl_id_iterator
   acl_id_iterator
   #(
      .WIDTH(32),
      .LOCAL_WIDTH_X(32),
      .LOCAL_WIDTH_Y(32),
      .LOCAL_WIDTH_Z(32),
      .ENABLE_TESSELLATION(1)
   )
   Krnl_GA_id_iter_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_GA_wg_disp_start_out),
      .valid_in(Krnl_GA_wg_disp_valid_out[0]),
      .stall_out(Krnl_GA_wg_disp_stall_in[0]),
      .stall_in(Krnl_GA_stall_out[0]),
      .valid_out(Krnl_GA_valid_in[0]),
      .group_id_in(Krnl_GA_wg_disp_group_id_out),
      .global_id_base_in(Krnl_GA_wg_disp_global_id_base_out),
      .local_size(Krnl_GA_local_size),
      .global_size(Krnl_GA_global_size),
      .local_id(Krnl_GA_local_id[0]),
      .global_id(Krnl_GA_global_id[0]),
      .group_id(Krnl_GA_group_id[0])
   );

   // INST Krnl_GA_inst_0 of Krnl_GA_top_wrapper_0
   Krnl_GA_top_wrapper_0 Krnl_GA_inst_0
   (
      .start(Krnl_GA_start_kernel_copy[0]),
      .kernel_arguments(Krnl_GA_kernel_arguments),
      .work_dim(Krnl_GA_work_dim),
      .global_offset(Krnl_GA_global_offset),
      .kernel_valid_out(Krnl_GA_valid_out[0]),
      .has_a_write_pending(Krnl_GA_pending_write[0]),
      .has_a_lsu_active(Krnl_GA_lsu_active[0]),
      .global_id(Krnl_GA_global_id[0]),
      .local_id(Krnl_GA_local_id[0]),
      .group_id(Krnl_GA_group_id[0]),
      .global_size(Krnl_GA_global_size),
      .local_size(Krnl_GA_local_size),
      .num_groups(Krnl_GA_num_groups),
      .workgroup_size(Krnl_GA_wg_size),
      .kernel_stall_out(Krnl_GA_stall_out[0]),
      .kernel_valid_in(Krnl_GA_valid_in[0]),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_st__inst0
      .avm_local_bb0_st__inst0_enable(avm_kernel_wr_enable[1]),
      .avm_local_bb0_st__inst0_read(avm_kernel_wr_read[1]),
      .avm_local_bb0_st__inst0_write(avm_kernel_wr_write[1]),
      .avm_local_bb0_st__inst0_burstcount(avm_kernel_wr_burstcount[1]),
      .avm_local_bb0_st__inst0_address(avm_kernel_wr_address[1]),
      .avm_local_bb0_st__inst0_writedata(avm_kernel_wr_writedata[1]),
      .avm_local_bb0_st__inst0_byteenable(avm_kernel_wr_byteenable[1]),
      .avm_local_bb0_st__inst0_waitrequest(avm_kernel_wr_waitrequest[1]),
      .avm_local_bb0_st__inst0_readdata(avm_kernel_wr_readdata[1]),
      .avm_local_bb0_st__inst0_readdatavalid(avm_kernel_wr_readdatavalid[1]),
      .avm_local_bb0_st__inst0_writeack(avm_kernel_wr_writeack[1]),
      // AVM avm_local_bb2_st__inst0
      .avm_local_bb2_st__inst0_enable(avm_kernel_wr_enable[2]),
      .avm_local_bb2_st__inst0_read(avm_kernel_wr_read[2]),
      .avm_local_bb2_st__inst0_write(avm_kernel_wr_write[2]),
      .avm_local_bb2_st__inst0_burstcount(avm_kernel_wr_burstcount[2]),
      .avm_local_bb2_st__inst0_address(avm_kernel_wr_address[2]),
      .avm_local_bb2_st__inst0_writedata(avm_kernel_wr_writedata[2]),
      .avm_local_bb2_st__inst0_byteenable(avm_kernel_wr_byteenable[2]),
      .avm_local_bb2_st__inst0_waitrequest(avm_kernel_wr_waitrequest[2]),
      .avm_local_bb2_st__inst0_readdata(avm_kernel_wr_readdata[2]),
      .avm_local_bb2_st__inst0_readdatavalid(avm_kernel_wr_readdatavalid[2]),
      .avm_local_bb2_st__inst0_writeack(avm_kernel_wr_writeack[2]),
      // AVM avm_local_bb2_st__u1_inst0
      .avm_local_bb2_st__u1_inst0_enable(avm_kernel_wr_enable[3]),
      .avm_local_bb2_st__u1_inst0_read(avm_kernel_wr_read[3]),
      .avm_local_bb2_st__u1_inst0_write(avm_kernel_wr_write[3]),
      .avm_local_bb2_st__u1_inst0_burstcount(avm_kernel_wr_burstcount[3]),
      .avm_local_bb2_st__u1_inst0_address(avm_kernel_wr_address[3]),
      .avm_local_bb2_st__u1_inst0_writedata(avm_kernel_wr_writedata[3]),
      .avm_local_bb2_st__u1_inst0_byteenable(avm_kernel_wr_byteenable[3]),
      .avm_local_bb2_st__u1_inst0_waitrequest(avm_kernel_wr_waitrequest[3]),
      .avm_local_bb2_st__u1_inst0_readdata(avm_kernel_wr_readdata[3]),
      .avm_local_bb2_st__u1_inst0_readdatavalid(avm_kernel_wr_readdatavalid[3]),
      .avm_local_bb2_st__u1_inst0_writeack(avm_kernel_wr_writeack[3]),
      // AVST avm_channel_id_chan_GA2Conf_genotype_write
      .avm_channel_id_chan_GA2Conf_genotype_write_valid(avm_channel_id_chan_GA2Conf_genotype_write_valid),
      .avm_channel_id_chan_GA2Conf_genotype_write_ready(avm_channel_id_chan_GA2Conf_genotype_write_ready),
      .avm_channel_id_chan_GA2Conf_genotype_write_data(avm_channel_id_chan_GA2Conf_genotype_write_data),
      .avm_channel_id_chan_GA2Conf_genotype_write_almostfull(avm_channel_id_chan_GA2Conf_genotype_write_almostfull),
      // AVST avm_channel_id_chan_GA2Conf_active_write
      .avm_channel_id_chan_GA2Conf_active_write_valid(avm_channel_id_chan_GA2Conf_active_write_valid),
      .avm_channel_id_chan_GA2Conf_active_write_ready(avm_channel_id_chan_GA2Conf_active_write_ready),
      .avm_channel_id_chan_GA2Conf_active_write_data(avm_channel_id_chan_GA2Conf_active_write_data),
      .avm_channel_id_chan_GA2Conf_active_write_almostfull(avm_channel_id_chan_GA2Conf_active_write_almostfull),
      // AVST avm_channel_id_chan_GA2Conf_cnt_write
      .avm_channel_id_chan_GA2Conf_cnt_write_valid(avm_channel_id_chan_GA2Conf_cnt_write_valid),
      .avm_channel_id_chan_GA2Conf_cnt_write_ready(avm_channel_id_chan_GA2Conf_cnt_write_ready),
      .avm_channel_id_chan_GA2Conf_cnt_write_data(avm_channel_id_chan_GA2Conf_cnt_write_data),
      .avm_channel_id_chan_GA2Conf_cnt_write_almostfull(avm_channel_id_chan_GA2Conf_cnt_write_almostfull),
      // AVST avm_channel_id_chan_GA2Conf_mode_write
      .avm_channel_id_chan_GA2Conf_mode_write_valid(avm_channel_id_chan_GA2Conf_mode_write_valid),
      .avm_channel_id_chan_GA2Conf_mode_write_ready(avm_channel_id_chan_GA2Conf_mode_write_ready),
      .avm_channel_id_chan_GA2Conf_mode_write_data(avm_channel_id_chan_GA2Conf_mode_write_data),
      .avm_channel_id_chan_GA2Conf_mode_write_almostfull(avm_channel_id_chan_GA2Conf_mode_write_almostfull),
      // AVM p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_GA_enable[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_GA_read[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_GA_write[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_GA_burstcount[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_GA_address[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_GA_writedata[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_GA_byteenable[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_GA_waitrequest[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_GA_readdata[0]),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_GA_readdatavalid[0])
   );

   // INST Krnl_GA_start_elem_inst_0 of acl_start_signal_chain_element
   acl_start_signal_chain_element Krnl_GA_start_elem_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start_in(Krnl_GA_start_chain[0]),
      .start_kernel(Krnl_GA_start_kernel_copy[0]),
      .start_finish_detector(Krnl_GA_start_task_fd[0]),
      .start_finish_chain_element(Krnl_GA_start_finish_element[0]),
      .start_chain()
   );

   assign Krnl_InterE_start_chain[0] = Krnl_InterE_start;
   assign Krnl_InterE_finish_chain[0] = 1'b1;
   assign Krnl_InterE_cra_pending_write = |Krnl_InterE_pending_write;
   assign Krnl_InterE_cra_lsu_active = |Krnl_InterE_lsu_active;
   assign Krnl_InterE_cra_valid_in = |Krnl_InterE_valid_in;
   assign Krnl_InterE_stall_in = 0;
   // INST Krnl_InterE_workgroup_dispatcher of acl_work_group_dispatcher
   acl_work_group_dispatcher
   #(
      .WIDTH(32),
      .NUM_COPIES(1),
      .RUN_FOREVER(0)
   )
   Krnl_InterE_workgroup_dispatcher
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_InterE_start),
      .num_groups(Krnl_InterE_num_groups),
      .local_size(Krnl_InterE_local_size),
      .stall_in(Krnl_InterE_wg_disp_stall_in),
      .valid_out(Krnl_InterE_wg_disp_valid_out),
      .group_id_out(Krnl_InterE_wg_disp_group_id_out),
      .global_id_base_out(Krnl_InterE_wg_disp_global_id_base_out),
      .start_out(Krnl_InterE_wg_disp_start_out),
      .dispatched_all_groups(Krnl_InterE_wg_disp_dispatched_all_groups)
   );

   // INST Krnl_InterE_finish_detector of acl_kernel_finish_detector
   acl_kernel_finish_detector
   #(
      .NUM_COPIES(1),
      .WG_SIZE_W(32),
      .GLOBAL_ID_W(32),
      .TESSELLATION_SIZE(19)
   )
   Krnl_InterE_finish_detector
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_InterE_start),
      .wg_size(Krnl_InterE_wg_size),
      .wg_dispatch_valid_out(Krnl_InterE_wg_disp_valid_out),
      .wg_dispatch_stall_in(Krnl_InterE_wg_disp_stall_in),
      .dispatched_all_groups(Krnl_InterE_wg_disp_dispatched_all_groups),
      .kernel_copy_valid_out(Krnl_InterE_valid_out),
      .kernel_copy_stall_in(Krnl_InterE_stall_in),
      .pending_writes(Krnl_InterE_cra_pending_write),
      .finish(Krnl_InterE_finish)
   );

   // INST Krnl_InterE_cra_slave_inst of Krnl_InterE_function_cra_slave
   Krnl_InterE_function_cra_slave Krnl_InterE_cra_slave_inst
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_InterE_start),
      .finish(Krnl_InterE_finish),
      .global_offset_0(Krnl_InterE_global_offset[0]),
      .global_offset_1(Krnl_InterE_global_offset[1]),
      .global_offset_2(Krnl_InterE_global_offset[2]),
      .work_dim(Krnl_InterE_work_dim),
      .has_a_lsu_active(Krnl_InterE_cra_lsu_active),
      .has_a_write_pending(Krnl_InterE_cra_pending_write),
      .valid_in(Krnl_InterE_cra_valid_in),
      .global_size_0(Krnl_InterE_global_size[0]),
      .global_size_1(Krnl_InterE_global_size[1]),
      .global_size_2(Krnl_InterE_global_size[2]),
      .num_groups_0(Krnl_InterE_num_groups[0]),
      .num_groups_1(Krnl_InterE_num_groups[1]),
      .num_groups_2(Krnl_InterE_num_groups[2]),
      .local_size_0(Krnl_InterE_local_size[0]),
      .local_size_1(Krnl_InterE_local_size[1]),
      .local_size_2(Krnl_InterE_local_size[2]),
      .workgroup_size(Krnl_InterE_wg_size),
      .kernel_arguments(Krnl_InterE_kernel_arguments),
      .cra_irq(kernel_irqs[2]),
      // AVS avs_cra
      .avs_cra_enable(avs_Krnl_InterE_cra_enable),
      .avs_cra_read(avs_Krnl_InterE_cra_read),
      .avs_cra_write(avs_Krnl_InterE_cra_write),
      .avs_cra_address(avs_Krnl_InterE_cra_address),
      .avs_cra_writedata(avs_Krnl_InterE_cra_writedata),
      .avs_cra_byteenable(avs_Krnl_InterE_cra_byteenable),
      .avs_cra_readdata(avs_Krnl_InterE_cra_readdata),
      .avs_cra_readdatavalid(avs_Krnl_InterE_cra_readdatavalid),
      .acl_counter_reset(counter_reset_Krnl_InterE),
      .acl_counter_init(counter_init_Krnl_InterE),
      .acl_counter_limit(counter_limit_Krnl_InterE),
      .acl_counter_size(counter_size_Krnl_InterE),
      .acl_counter_full(counter_full_Krnl_InterE)
   );

   // INST Krnl_InterE_id_iter_inst_0 of acl_id_iterator
   acl_id_iterator
   #(
      .WIDTH(32),
      .LOCAL_WIDTH_X(1),
      .LOCAL_WIDTH_Y(1),
      .LOCAL_WIDTH_Z(1),
      .ENABLE_TESSELLATION(1)
   )
   Krnl_InterE_id_iter_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_InterE_wg_disp_start_out),
      .valid_in(Krnl_InterE_wg_disp_valid_out[0]),
      .stall_out(Krnl_InterE_wg_disp_stall_in[0]),
      .stall_in(Krnl_InterE_stall_out[0]),
      .valid_out(Krnl_InterE_valid_in[0]),
      .group_id_in(Krnl_InterE_wg_disp_group_id_out),
      .global_id_base_in(Krnl_InterE_wg_disp_global_id_base_out),
      .local_size(Krnl_InterE_local_size),
      .global_size(Krnl_InterE_global_size),
      .local_id(Krnl_InterE_local_id[0]),
      .global_id(Krnl_InterE_global_id[0]),
      .group_id(Krnl_InterE_group_id[0])
   );

   // INST Krnl_InterE_inst_0 of Krnl_InterE_top_wrapper_0
   Krnl_InterE_top_wrapper_0 Krnl_InterE_inst_0
   (
      .start(Krnl_InterE_start_kernel_copy[0]),
      .kernel_arguments(Krnl_InterE_kernel_arguments),
      .work_dim(Krnl_InterE_work_dim),
      .global_offset(Krnl_InterE_global_offset),
      .kernel_valid_out(Krnl_InterE_valid_out[0]),
      .has_a_write_pending(Krnl_InterE_pending_write[0]),
      .has_a_lsu_active(Krnl_InterE_lsu_active[0]),
      .global_id(Krnl_InterE_global_id[0]),
      .local_id(Krnl_InterE_local_id[0]),
      .group_id(Krnl_InterE_group_id[0]),
      .global_size(Krnl_InterE_global_size),
      .local_size(Krnl_InterE_local_size),
      .num_groups(Krnl_InterE_num_groups),
      .workgroup_size(Krnl_InterE_wg_size),
      .kernel_stall_out(Krnl_InterE_stall_out[0]),
      .kernel_valid_in(Krnl_InterE_valid_in[0]),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable(avm_kernel_rd_enable[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read(avm_kernel_rd_read[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write(avm_kernel_rd_write[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount(avm_kernel_rd_burstcount[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address(avm_kernel_rd_address[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata(avm_kernel_rd_writedata[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable(avm_kernel_rd_byteenable[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest(avm_kernel_rd_waitrequest[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata(avm_kernel_rd_readdata[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid(avm_kernel_rd_readdatavalid[8]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack(avm_kernel_rd_writeack[8]),
      // AVM avm_local_bb3_st__inst0
      .avm_local_bb3_st__inst0_enable(avm_kernel_wr_enable[4]),
      .avm_local_bb3_st__inst0_read(avm_kernel_wr_read[4]),
      .avm_local_bb3_st__inst0_write(avm_kernel_wr_write[4]),
      .avm_local_bb3_st__inst0_burstcount(avm_kernel_wr_burstcount[4]),
      .avm_local_bb3_st__inst0_address(avm_kernel_wr_address[4]),
      .avm_local_bb3_st__inst0_writedata(avm_kernel_wr_writedata[4]),
      .avm_local_bb3_st__inst0_byteenable(avm_kernel_wr_byteenable[4]),
      .avm_local_bb3_st__inst0_waitrequest(avm_kernel_wr_waitrequest[4]),
      .avm_local_bb3_st__inst0_readdata(avm_kernel_wr_readdata[4]),
      .avm_local_bb3_st__inst0_readdatavalid(avm_kernel_wr_readdatavalid[4]),
      .avm_local_bb3_st__inst0_writeack(avm_kernel_wr_writeack[4]),
      // AVM avm_local_bb4_ld__inst0
      .avm_local_bb4_ld__inst0_enable(avm_kernel_rd_enable[9]),
      .avm_local_bb4_ld__inst0_read(avm_kernel_rd_read[9]),
      .avm_local_bb4_ld__inst0_write(avm_kernel_rd_write[9]),
      .avm_local_bb4_ld__inst0_burstcount(avm_kernel_rd_burstcount[9]),
      .avm_local_bb4_ld__inst0_address(avm_kernel_rd_address[9]),
      .avm_local_bb4_ld__inst0_writedata(avm_kernel_rd_writedata[9]),
      .avm_local_bb4_ld__inst0_byteenable(avm_kernel_rd_byteenable[9]),
      .avm_local_bb4_ld__inst0_waitrequest(avm_kernel_rd_waitrequest[9]),
      .avm_local_bb4_ld__inst0_readdata(avm_kernel_rd_readdata[9]),
      .avm_local_bb4_ld__inst0_readdatavalid(avm_kernel_rd_readdatavalid[9]),
      .avm_local_bb4_ld__inst0_writeack(avm_kernel_rd_writeack[9]),
      // AVM avm_local_bb4_ld__u11_inst0
      .avm_local_bb4_ld__u11_inst0_enable(avm_kernel_rd_enable[10]),
      .avm_local_bb4_ld__u11_inst0_read(avm_kernel_rd_read[10]),
      .avm_local_bb4_ld__u11_inst0_write(avm_kernel_rd_write[10]),
      .avm_local_bb4_ld__u11_inst0_burstcount(avm_kernel_rd_burstcount[10]),
      .avm_local_bb4_ld__u11_inst0_address(avm_kernel_rd_address[10]),
      .avm_local_bb4_ld__u11_inst0_writedata(avm_kernel_rd_writedata[10]),
      .avm_local_bb4_ld__u11_inst0_byteenable(avm_kernel_rd_byteenable[10]),
      .avm_local_bb4_ld__u11_inst0_waitrequest(avm_kernel_rd_waitrequest[10]),
      .avm_local_bb4_ld__u11_inst0_readdata(avm_kernel_rd_readdata[10]),
      .avm_local_bb4_ld__u11_inst0_readdatavalid(avm_kernel_rd_readdatavalid[10]),
      .avm_local_bb4_ld__u11_inst0_writeack(avm_kernel_rd_writeack[10]),
      // AVM avm_local_bb4_ld__u12_inst0
      .avm_local_bb4_ld__u12_inst0_enable(avm_kernel_rd_enable[11]),
      .avm_local_bb4_ld__u12_inst0_read(avm_kernel_rd_read[11]),
      .avm_local_bb4_ld__u12_inst0_write(avm_kernel_rd_write[11]),
      .avm_local_bb4_ld__u12_inst0_burstcount(avm_kernel_rd_burstcount[11]),
      .avm_local_bb4_ld__u12_inst0_address(avm_kernel_rd_address[11]),
      .avm_local_bb4_ld__u12_inst0_writedata(avm_kernel_rd_writedata[11]),
      .avm_local_bb4_ld__u12_inst0_byteenable(avm_kernel_rd_byteenable[11]),
      .avm_local_bb4_ld__u12_inst0_waitrequest(avm_kernel_rd_waitrequest[11]),
      .avm_local_bb4_ld__u12_inst0_readdata(avm_kernel_rd_readdata[11]),
      .avm_local_bb4_ld__u12_inst0_readdatavalid(avm_kernel_rd_readdatavalid[11]),
      .avm_local_bb4_ld__u12_inst0_writeack(avm_kernel_rd_writeack[11]),
      // AVM avm_local_bb4_ld__u13_inst0
      .avm_local_bb4_ld__u13_inst0_enable(avm_kernel_rd_enable[12]),
      .avm_local_bb4_ld__u13_inst0_read(avm_kernel_rd_read[12]),
      .avm_local_bb4_ld__u13_inst0_write(avm_kernel_rd_write[12]),
      .avm_local_bb4_ld__u13_inst0_burstcount(avm_kernel_rd_burstcount[12]),
      .avm_local_bb4_ld__u13_inst0_address(avm_kernel_rd_address[12]),
      .avm_local_bb4_ld__u13_inst0_writedata(avm_kernel_rd_writedata[12]),
      .avm_local_bb4_ld__u13_inst0_byteenable(avm_kernel_rd_byteenable[12]),
      .avm_local_bb4_ld__u13_inst0_waitrequest(avm_kernel_rd_waitrequest[12]),
      .avm_local_bb4_ld__u13_inst0_readdata(avm_kernel_rd_readdata[12]),
      .avm_local_bb4_ld__u13_inst0_readdatavalid(avm_kernel_rd_readdatavalid[12]),
      .avm_local_bb4_ld__u13_inst0_writeack(avm_kernel_rd_writeack[12]),
      // AVM avm_local_bb4_ld__u14_inst0
      .avm_local_bb4_ld__u14_inst0_enable(avm_kernel_rd_enable[13]),
      .avm_local_bb4_ld__u14_inst0_read(avm_kernel_rd_read[13]),
      .avm_local_bb4_ld__u14_inst0_write(avm_kernel_rd_write[13]),
      .avm_local_bb4_ld__u14_inst0_burstcount(avm_kernel_rd_burstcount[13]),
      .avm_local_bb4_ld__u14_inst0_address(avm_kernel_rd_address[13]),
      .avm_local_bb4_ld__u14_inst0_writedata(avm_kernel_rd_writedata[13]),
      .avm_local_bb4_ld__u14_inst0_byteenable(avm_kernel_rd_byteenable[13]),
      .avm_local_bb4_ld__u14_inst0_waitrequest(avm_kernel_rd_waitrequest[13]),
      .avm_local_bb4_ld__u14_inst0_readdata(avm_kernel_rd_readdata[13]),
      .avm_local_bb4_ld__u14_inst0_readdatavalid(avm_kernel_rd_readdatavalid[13]),
      .avm_local_bb4_ld__u14_inst0_writeack(avm_kernel_rd_writeack[13]),
      // AVM avm_local_bb4_ld__u15_inst0
      .avm_local_bb4_ld__u15_inst0_enable(avm_kernel_rd_enable[14]),
      .avm_local_bb4_ld__u15_inst0_read(avm_kernel_rd_read[14]),
      .avm_local_bb4_ld__u15_inst0_write(avm_kernel_rd_write[14]),
      .avm_local_bb4_ld__u15_inst0_burstcount(avm_kernel_rd_burstcount[14]),
      .avm_local_bb4_ld__u15_inst0_address(avm_kernel_rd_address[14]),
      .avm_local_bb4_ld__u15_inst0_writedata(avm_kernel_rd_writedata[14]),
      .avm_local_bb4_ld__u15_inst0_byteenable(avm_kernel_rd_byteenable[14]),
      .avm_local_bb4_ld__u15_inst0_waitrequest(avm_kernel_rd_waitrequest[14]),
      .avm_local_bb4_ld__u15_inst0_readdata(avm_kernel_rd_readdata[14]),
      .avm_local_bb4_ld__u15_inst0_readdatavalid(avm_kernel_rd_readdatavalid[14]),
      .avm_local_bb4_ld__u15_inst0_writeack(avm_kernel_rd_writeack[14]),
      // AVM avm_local_bb4_ld__u16_inst0
      .avm_local_bb4_ld__u16_inst0_enable(avm_kernel_rd_enable[15]),
      .avm_local_bb4_ld__u16_inst0_read(avm_kernel_rd_read[15]),
      .avm_local_bb4_ld__u16_inst0_write(avm_kernel_rd_write[15]),
      .avm_local_bb4_ld__u16_inst0_burstcount(avm_kernel_rd_burstcount[15]),
      .avm_local_bb4_ld__u16_inst0_address(avm_kernel_rd_address[15]),
      .avm_local_bb4_ld__u16_inst0_writedata(avm_kernel_rd_writedata[15]),
      .avm_local_bb4_ld__u16_inst0_byteenable(avm_kernel_rd_byteenable[15]),
      .avm_local_bb4_ld__u16_inst0_waitrequest(avm_kernel_rd_waitrequest[15]),
      .avm_local_bb4_ld__u16_inst0_readdata(avm_kernel_rd_readdata[15]),
      .avm_local_bb4_ld__u16_inst0_readdatavalid(avm_kernel_rd_readdatavalid[15]),
      .avm_local_bb4_ld__u16_inst0_writeack(avm_kernel_rd_writeack[15]),
      // AVM avm_local_bb4_ld__u17_inst0
      .avm_local_bb4_ld__u17_inst0_enable(avm_kernel_rd_enable[16]),
      .avm_local_bb4_ld__u17_inst0_read(avm_kernel_rd_read[16]),
      .avm_local_bb4_ld__u17_inst0_write(avm_kernel_rd_write[16]),
      .avm_local_bb4_ld__u17_inst0_burstcount(avm_kernel_rd_burstcount[16]),
      .avm_local_bb4_ld__u17_inst0_address(avm_kernel_rd_address[16]),
      .avm_local_bb4_ld__u17_inst0_writedata(avm_kernel_rd_writedata[16]),
      .avm_local_bb4_ld__u17_inst0_byteenable(avm_kernel_rd_byteenable[16]),
      .avm_local_bb4_ld__u17_inst0_waitrequest(avm_kernel_rd_waitrequest[16]),
      .avm_local_bb4_ld__u17_inst0_readdata(avm_kernel_rd_readdata[16]),
      .avm_local_bb4_ld__u17_inst0_readdatavalid(avm_kernel_rd_readdatavalid[16]),
      .avm_local_bb4_ld__u17_inst0_writeack(avm_kernel_rd_writeack[16]),
      // AVM avm_local_bb4_ld__u18_inst0
      .avm_local_bb4_ld__u18_inst0_enable(avm_kernel_rd_enable[17]),
      .avm_local_bb4_ld__u18_inst0_read(avm_kernel_rd_read[17]),
      .avm_local_bb4_ld__u18_inst0_write(avm_kernel_rd_write[17]),
      .avm_local_bb4_ld__u18_inst0_burstcount(avm_kernel_rd_burstcount[17]),
      .avm_local_bb4_ld__u18_inst0_address(avm_kernel_rd_address[17]),
      .avm_local_bb4_ld__u18_inst0_writedata(avm_kernel_rd_writedata[17]),
      .avm_local_bb4_ld__u18_inst0_byteenable(avm_kernel_rd_byteenable[17]),
      .avm_local_bb4_ld__u18_inst0_waitrequest(avm_kernel_rd_waitrequest[17]),
      .avm_local_bb4_ld__u18_inst0_readdata(avm_kernel_rd_readdata[17]),
      .avm_local_bb4_ld__u18_inst0_readdatavalid(avm_kernel_rd_readdatavalid[17]),
      .avm_local_bb4_ld__u18_inst0_writeack(avm_kernel_rd_writeack[17]),
      // AVM avm_local_bb4_ld__u19_inst0
      .avm_local_bb4_ld__u19_inst0_enable(avm_kernel_rd_enable[18]),
      .avm_local_bb4_ld__u19_inst0_read(avm_kernel_rd_read[18]),
      .avm_local_bb4_ld__u19_inst0_write(avm_kernel_rd_write[18]),
      .avm_local_bb4_ld__u19_inst0_burstcount(avm_kernel_rd_burstcount[18]),
      .avm_local_bb4_ld__u19_inst0_address(avm_kernel_rd_address[18]),
      .avm_local_bb4_ld__u19_inst0_writedata(avm_kernel_rd_writedata[18]),
      .avm_local_bb4_ld__u19_inst0_byteenable(avm_kernel_rd_byteenable[18]),
      .avm_local_bb4_ld__u19_inst0_waitrequest(avm_kernel_rd_waitrequest[18]),
      .avm_local_bb4_ld__u19_inst0_readdata(avm_kernel_rd_readdata[18]),
      .avm_local_bb4_ld__u19_inst0_readdatavalid(avm_kernel_rd_readdatavalid[18]),
      .avm_local_bb4_ld__u19_inst0_writeack(avm_kernel_rd_writeack[18]),
      // AVM avm_local_bb4_ld__u20_inst0
      .avm_local_bb4_ld__u20_inst0_enable(avm_kernel_rd_enable[19]),
      .avm_local_bb4_ld__u20_inst0_read(avm_kernel_rd_read[19]),
      .avm_local_bb4_ld__u20_inst0_write(avm_kernel_rd_write[19]),
      .avm_local_bb4_ld__u20_inst0_burstcount(avm_kernel_rd_burstcount[19]),
      .avm_local_bb4_ld__u20_inst0_address(avm_kernel_rd_address[19]),
      .avm_local_bb4_ld__u20_inst0_writedata(avm_kernel_rd_writedata[19]),
      .avm_local_bb4_ld__u20_inst0_byteenable(avm_kernel_rd_byteenable[19]),
      .avm_local_bb4_ld__u20_inst0_waitrequest(avm_kernel_rd_waitrequest[19]),
      .avm_local_bb4_ld__u20_inst0_readdata(avm_kernel_rd_readdata[19]),
      .avm_local_bb4_ld__u20_inst0_readdatavalid(avm_kernel_rd_readdatavalid[19]),
      .avm_local_bb4_ld__u20_inst0_writeack(avm_kernel_rd_writeack[19]),
      // AVM avm_local_bb4_ld__u21_inst0
      .avm_local_bb4_ld__u21_inst0_enable(avm_kernel_rd_enable[20]),
      .avm_local_bb4_ld__u21_inst0_read(avm_kernel_rd_read[20]),
      .avm_local_bb4_ld__u21_inst0_write(avm_kernel_rd_write[20]),
      .avm_local_bb4_ld__u21_inst0_burstcount(avm_kernel_rd_burstcount[20]),
      .avm_local_bb4_ld__u21_inst0_address(avm_kernel_rd_address[20]),
      .avm_local_bb4_ld__u21_inst0_writedata(avm_kernel_rd_writedata[20]),
      .avm_local_bb4_ld__u21_inst0_byteenable(avm_kernel_rd_byteenable[20]),
      .avm_local_bb4_ld__u21_inst0_waitrequest(avm_kernel_rd_waitrequest[20]),
      .avm_local_bb4_ld__u21_inst0_readdata(avm_kernel_rd_readdata[20]),
      .avm_local_bb4_ld__u21_inst0_readdatavalid(avm_kernel_rd_readdatavalid[20]),
      .avm_local_bb4_ld__u21_inst0_writeack(avm_kernel_rd_writeack[20]),
      // AVM avm_local_bb4_ld__u22_inst0
      .avm_local_bb4_ld__u22_inst0_enable(avm_kernel_rd_enable[21]),
      .avm_local_bb4_ld__u22_inst0_read(avm_kernel_rd_read[21]),
      .avm_local_bb4_ld__u22_inst0_write(avm_kernel_rd_write[21]),
      .avm_local_bb4_ld__u22_inst0_burstcount(avm_kernel_rd_burstcount[21]),
      .avm_local_bb4_ld__u22_inst0_address(avm_kernel_rd_address[21]),
      .avm_local_bb4_ld__u22_inst0_writedata(avm_kernel_rd_writedata[21]),
      .avm_local_bb4_ld__u22_inst0_byteenable(avm_kernel_rd_byteenable[21]),
      .avm_local_bb4_ld__u22_inst0_waitrequest(avm_kernel_rd_waitrequest[21]),
      .avm_local_bb4_ld__u22_inst0_readdata(avm_kernel_rd_readdata[21]),
      .avm_local_bb4_ld__u22_inst0_readdatavalid(avm_kernel_rd_readdatavalid[21]),
      .avm_local_bb4_ld__u22_inst0_writeack(avm_kernel_rd_writeack[21]),
      // AVM avm_local_bb4_ld__u23_inst0
      .avm_local_bb4_ld__u23_inst0_enable(avm_kernel_rd_enable[22]),
      .avm_local_bb4_ld__u23_inst0_read(avm_kernel_rd_read[22]),
      .avm_local_bb4_ld__u23_inst0_write(avm_kernel_rd_write[22]),
      .avm_local_bb4_ld__u23_inst0_burstcount(avm_kernel_rd_burstcount[22]),
      .avm_local_bb4_ld__u23_inst0_address(avm_kernel_rd_address[22]),
      .avm_local_bb4_ld__u23_inst0_writedata(avm_kernel_rd_writedata[22]),
      .avm_local_bb4_ld__u23_inst0_byteenable(avm_kernel_rd_byteenable[22]),
      .avm_local_bb4_ld__u23_inst0_waitrequest(avm_kernel_rd_waitrequest[22]),
      .avm_local_bb4_ld__u23_inst0_readdata(avm_kernel_rd_readdata[22]),
      .avm_local_bb4_ld__u23_inst0_readdatavalid(avm_kernel_rd_readdatavalid[22]),
      .avm_local_bb4_ld__u23_inst0_writeack(avm_kernel_rd_writeack[22]),
      // AVM avm_local_bb4_ld__u24_inst0
      .avm_local_bb4_ld__u24_inst0_enable(avm_kernel_rd_enable[23]),
      .avm_local_bb4_ld__u24_inst0_read(avm_kernel_rd_read[23]),
      .avm_local_bb4_ld__u24_inst0_write(avm_kernel_rd_write[23]),
      .avm_local_bb4_ld__u24_inst0_burstcount(avm_kernel_rd_burstcount[23]),
      .avm_local_bb4_ld__u24_inst0_address(avm_kernel_rd_address[23]),
      .avm_local_bb4_ld__u24_inst0_writedata(avm_kernel_rd_writedata[23]),
      .avm_local_bb4_ld__u24_inst0_byteenable(avm_kernel_rd_byteenable[23]),
      .avm_local_bb4_ld__u24_inst0_waitrequest(avm_kernel_rd_waitrequest[23]),
      .avm_local_bb4_ld__u24_inst0_readdata(avm_kernel_rd_readdata[23]),
      .avm_local_bb4_ld__u24_inst0_readdatavalid(avm_kernel_rd_readdatavalid[23]),
      .avm_local_bb4_ld__u24_inst0_writeack(avm_kernel_rd_writeack[23]),
      // AVM avm_local_bb4_ld__u25_inst0
      .avm_local_bb4_ld__u25_inst0_enable(avm_kernel_rd_enable[24]),
      .avm_local_bb4_ld__u25_inst0_read(avm_kernel_rd_read[24]),
      .avm_local_bb4_ld__u25_inst0_write(avm_kernel_rd_write[24]),
      .avm_local_bb4_ld__u25_inst0_burstcount(avm_kernel_rd_burstcount[24]),
      .avm_local_bb4_ld__u25_inst0_address(avm_kernel_rd_address[24]),
      .avm_local_bb4_ld__u25_inst0_writedata(avm_kernel_rd_writedata[24]),
      .avm_local_bb4_ld__u25_inst0_byteenable(avm_kernel_rd_byteenable[24]),
      .avm_local_bb4_ld__u25_inst0_waitrequest(avm_kernel_rd_waitrequest[24]),
      .avm_local_bb4_ld__u25_inst0_readdata(avm_kernel_rd_readdata[24]),
      .avm_local_bb4_ld__u25_inst0_readdatavalid(avm_kernel_rd_readdatavalid[24]),
      .avm_local_bb4_ld__u25_inst0_writeack(avm_kernel_rd_writeack[24]),
      // AVM avm_local_bb4_ld__u26_inst0
      .avm_local_bb4_ld__u26_inst0_enable(avm_kernel_rd_enable[25]),
      .avm_local_bb4_ld__u26_inst0_read(avm_kernel_rd_read[25]),
      .avm_local_bb4_ld__u26_inst0_write(avm_kernel_rd_write[25]),
      .avm_local_bb4_ld__u26_inst0_burstcount(avm_kernel_rd_burstcount[25]),
      .avm_local_bb4_ld__u26_inst0_address(avm_kernel_rd_address[25]),
      .avm_local_bb4_ld__u26_inst0_writedata(avm_kernel_rd_writedata[25]),
      .avm_local_bb4_ld__u26_inst0_byteenable(avm_kernel_rd_byteenable[25]),
      .avm_local_bb4_ld__u26_inst0_waitrequest(avm_kernel_rd_waitrequest[25]),
      .avm_local_bb4_ld__u26_inst0_readdata(avm_kernel_rd_readdata[25]),
      .avm_local_bb4_ld__u26_inst0_readdatavalid(avm_kernel_rd_readdatavalid[25]),
      .avm_local_bb4_ld__u26_inst0_writeack(avm_kernel_rd_writeack[25]),
      // AVM avm_local_bb4_ld__u27_inst0
      .avm_local_bb4_ld__u27_inst0_enable(avm_kernel_rd_enable[26]),
      .avm_local_bb4_ld__u27_inst0_read(avm_kernel_rd_read[26]),
      .avm_local_bb4_ld__u27_inst0_write(avm_kernel_rd_write[26]),
      .avm_local_bb4_ld__u27_inst0_burstcount(avm_kernel_rd_burstcount[26]),
      .avm_local_bb4_ld__u27_inst0_address(avm_kernel_rd_address[26]),
      .avm_local_bb4_ld__u27_inst0_writedata(avm_kernel_rd_writedata[26]),
      .avm_local_bb4_ld__u27_inst0_byteenable(avm_kernel_rd_byteenable[26]),
      .avm_local_bb4_ld__u27_inst0_waitrequest(avm_kernel_rd_waitrequest[26]),
      .avm_local_bb4_ld__u27_inst0_readdata(avm_kernel_rd_readdata[26]),
      .avm_local_bb4_ld__u27_inst0_readdatavalid(avm_kernel_rd_readdatavalid[26]),
      .avm_local_bb4_ld__u27_inst0_writeack(avm_kernel_rd_writeack[26]),
      // AVM avm_local_bb4_ld__u28_inst0
      .avm_local_bb4_ld__u28_inst0_enable(avm_kernel_rd_enable[27]),
      .avm_local_bb4_ld__u28_inst0_read(avm_kernel_rd_read[27]),
      .avm_local_bb4_ld__u28_inst0_write(avm_kernel_rd_write[27]),
      .avm_local_bb4_ld__u28_inst0_burstcount(avm_kernel_rd_burstcount[27]),
      .avm_local_bb4_ld__u28_inst0_address(avm_kernel_rd_address[27]),
      .avm_local_bb4_ld__u28_inst0_writedata(avm_kernel_rd_writedata[27]),
      .avm_local_bb4_ld__u28_inst0_byteenable(avm_kernel_rd_byteenable[27]),
      .avm_local_bb4_ld__u28_inst0_waitrequest(avm_kernel_rd_waitrequest[27]),
      .avm_local_bb4_ld__u28_inst0_readdata(avm_kernel_rd_readdata[27]),
      .avm_local_bb4_ld__u28_inst0_readdatavalid(avm_kernel_rd_readdatavalid[27]),
      .avm_local_bb4_ld__u28_inst0_writeack(avm_kernel_rd_writeack[27]),
      // AVM avm_local_bb4_ld__u29_inst0
      .avm_local_bb4_ld__u29_inst0_enable(avm_kernel_rd_enable[28]),
      .avm_local_bb4_ld__u29_inst0_read(avm_kernel_rd_read[28]),
      .avm_local_bb4_ld__u29_inst0_write(avm_kernel_rd_write[28]),
      .avm_local_bb4_ld__u29_inst0_burstcount(avm_kernel_rd_burstcount[28]),
      .avm_local_bb4_ld__u29_inst0_address(avm_kernel_rd_address[28]),
      .avm_local_bb4_ld__u29_inst0_writedata(avm_kernel_rd_writedata[28]),
      .avm_local_bb4_ld__u29_inst0_byteenable(avm_kernel_rd_byteenable[28]),
      .avm_local_bb4_ld__u29_inst0_waitrequest(avm_kernel_rd_waitrequest[28]),
      .avm_local_bb4_ld__u29_inst0_readdata(avm_kernel_rd_readdata[28]),
      .avm_local_bb4_ld__u29_inst0_readdatavalid(avm_kernel_rd_readdatavalid[28]),
      .avm_local_bb4_ld__u29_inst0_writeack(avm_kernel_rd_writeack[28]),
      // AVM avm_local_bb4_ld__u30_inst0
      .avm_local_bb4_ld__u30_inst0_enable(avm_kernel_rd_enable[29]),
      .avm_local_bb4_ld__u30_inst0_read(avm_kernel_rd_read[29]),
      .avm_local_bb4_ld__u30_inst0_write(avm_kernel_rd_write[29]),
      .avm_local_bb4_ld__u30_inst0_burstcount(avm_kernel_rd_burstcount[29]),
      .avm_local_bb4_ld__u30_inst0_address(avm_kernel_rd_address[29]),
      .avm_local_bb4_ld__u30_inst0_writedata(avm_kernel_rd_writedata[29]),
      .avm_local_bb4_ld__u30_inst0_byteenable(avm_kernel_rd_byteenable[29]),
      .avm_local_bb4_ld__u30_inst0_waitrequest(avm_kernel_rd_waitrequest[29]),
      .avm_local_bb4_ld__u30_inst0_readdata(avm_kernel_rd_readdata[29]),
      .avm_local_bb4_ld__u30_inst0_readdatavalid(avm_kernel_rd_readdatavalid[29]),
      .avm_local_bb4_ld__u30_inst0_writeack(avm_kernel_rd_writeack[29]),
      // AVM avm_local_bb4_ld__u31_inst0
      .avm_local_bb4_ld__u31_inst0_enable(avm_kernel_rd_enable[30]),
      .avm_local_bb4_ld__u31_inst0_read(avm_kernel_rd_read[30]),
      .avm_local_bb4_ld__u31_inst0_write(avm_kernel_rd_write[30]),
      .avm_local_bb4_ld__u31_inst0_burstcount(avm_kernel_rd_burstcount[30]),
      .avm_local_bb4_ld__u31_inst0_address(avm_kernel_rd_address[30]),
      .avm_local_bb4_ld__u31_inst0_writedata(avm_kernel_rd_writedata[30]),
      .avm_local_bb4_ld__u31_inst0_byteenable(avm_kernel_rd_byteenable[30]),
      .avm_local_bb4_ld__u31_inst0_waitrequest(avm_kernel_rd_waitrequest[30]),
      .avm_local_bb4_ld__u31_inst0_readdata(avm_kernel_rd_readdata[30]),
      .avm_local_bb4_ld__u31_inst0_readdatavalid(avm_kernel_rd_readdatavalid[30]),
      .avm_local_bb4_ld__u31_inst0_writeack(avm_kernel_rd_writeack[30]),
      // AVM avm_local_bb4_ld__u32_inst0
      .avm_local_bb4_ld__u32_inst0_enable(avm_kernel_rd_enable[31]),
      .avm_local_bb4_ld__u32_inst0_read(avm_kernel_rd_read[31]),
      .avm_local_bb4_ld__u32_inst0_write(avm_kernel_rd_write[31]),
      .avm_local_bb4_ld__u32_inst0_burstcount(avm_kernel_rd_burstcount[31]),
      .avm_local_bb4_ld__u32_inst0_address(avm_kernel_rd_address[31]),
      .avm_local_bb4_ld__u32_inst0_writedata(avm_kernel_rd_writedata[31]),
      .avm_local_bb4_ld__u32_inst0_byteenable(avm_kernel_rd_byteenable[31]),
      .avm_local_bb4_ld__u32_inst0_waitrequest(avm_kernel_rd_waitrequest[31]),
      .avm_local_bb4_ld__u32_inst0_readdata(avm_kernel_rd_readdata[31]),
      .avm_local_bb4_ld__u32_inst0_readdatavalid(avm_kernel_rd_readdatavalid[31]),
      .avm_local_bb4_ld__u32_inst0_writeack(avm_kernel_rd_writeack[31]),
      // AVM avm_local_bb4_ld__u33_inst0
      .avm_local_bb4_ld__u33_inst0_enable(avm_kernel_rd_enable[32]),
      .avm_local_bb4_ld__u33_inst0_read(avm_kernel_rd_read[32]),
      .avm_local_bb4_ld__u33_inst0_write(avm_kernel_rd_write[32]),
      .avm_local_bb4_ld__u33_inst0_burstcount(avm_kernel_rd_burstcount[32]),
      .avm_local_bb4_ld__u33_inst0_address(avm_kernel_rd_address[32]),
      .avm_local_bb4_ld__u33_inst0_writedata(avm_kernel_rd_writedata[32]),
      .avm_local_bb4_ld__u33_inst0_byteenable(avm_kernel_rd_byteenable[32]),
      .avm_local_bb4_ld__u33_inst0_waitrequest(avm_kernel_rd_waitrequest[32]),
      .avm_local_bb4_ld__u33_inst0_readdata(avm_kernel_rd_readdata[32]),
      .avm_local_bb4_ld__u33_inst0_readdatavalid(avm_kernel_rd_readdatavalid[32]),
      .avm_local_bb4_ld__u33_inst0_writeack(avm_kernel_rd_writeack[32]),
      // AVM avm_local_bb4_ld__u34_inst0
      .avm_local_bb4_ld__u34_inst0_enable(avm_kernel_rd_enable[33]),
      .avm_local_bb4_ld__u34_inst0_read(avm_kernel_rd_read[33]),
      .avm_local_bb4_ld__u34_inst0_write(avm_kernel_rd_write[33]),
      .avm_local_bb4_ld__u34_inst0_burstcount(avm_kernel_rd_burstcount[33]),
      .avm_local_bb4_ld__u34_inst0_address(avm_kernel_rd_address[33]),
      .avm_local_bb4_ld__u34_inst0_writedata(avm_kernel_rd_writedata[33]),
      .avm_local_bb4_ld__u34_inst0_byteenable(avm_kernel_rd_byteenable[33]),
      .avm_local_bb4_ld__u34_inst0_waitrequest(avm_kernel_rd_waitrequest[33]),
      .avm_local_bb4_ld__u34_inst0_readdata(avm_kernel_rd_readdata[33]),
      .avm_local_bb4_ld__u34_inst0_readdatavalid(avm_kernel_rd_readdatavalid[33]),
      .avm_local_bb4_ld__u34_inst0_writeack(avm_kernel_rd_writeack[33]),
      // AVM avm_local_bb4_ld__u7_inst0
      .avm_local_bb4_ld__u7_inst0_enable(avm_kernel_rd_enable[34]),
      .avm_local_bb4_ld__u7_inst0_read(avm_kernel_rd_read[34]),
      .avm_local_bb4_ld__u7_inst0_write(avm_kernel_rd_write[34]),
      .avm_local_bb4_ld__u7_inst0_burstcount(avm_kernel_rd_burstcount[34]),
      .avm_local_bb4_ld__u7_inst0_address(avm_kernel_rd_address[34]),
      .avm_local_bb4_ld__u7_inst0_writedata(avm_kernel_rd_writedata[34]),
      .avm_local_bb4_ld__u7_inst0_byteenable(avm_kernel_rd_byteenable[34]),
      .avm_local_bb4_ld__u7_inst0_waitrequest(avm_kernel_rd_waitrequest[34]),
      .avm_local_bb4_ld__u7_inst0_readdata(avm_kernel_rd_readdata[34]),
      .avm_local_bb4_ld__u7_inst0_readdatavalid(avm_kernel_rd_readdatavalid[34]),
      .avm_local_bb4_ld__u7_inst0_writeack(avm_kernel_rd_writeack[34]),
      // AVST avm_channel_id_chan_Conf2Intere_active_read
      .avm_channel_id_chan_Conf2Intere_active_read_valid(avm_channel_id_chan_Conf2Intere_active_read_valid),
      .avm_channel_id_chan_Conf2Intere_active_read_ready(avm_channel_id_chan_Conf2Intere_active_read_ready),
      .avm_channel_id_chan_Conf2Intere_active_read_data(avm_channel_id_chan_Conf2Intere_active_read_data),
      // AVST avm_channel_id_chan_Conf2Intere_cnt_read
      .avm_channel_id_chan_Conf2Intere_cnt_read_valid(avm_channel_id_chan_Conf2Intere_cnt_read_valid),
      .avm_channel_id_chan_Conf2Intere_cnt_read_ready(avm_channel_id_chan_Conf2Intere_cnt_read_ready),
      .avm_channel_id_chan_Conf2Intere_cnt_read_data(avm_channel_id_chan_Conf2Intere_cnt_read_data),
      // AVST avm_channel_id_chan_Conf2Intere_mode_read
      .avm_channel_id_chan_Conf2Intere_mode_read_valid(avm_channel_id_chan_Conf2Intere_mode_read_valid),
      .avm_channel_id_chan_Conf2Intere_mode_read_ready(avm_channel_id_chan_Conf2Intere_mode_read_ready),
      .avm_channel_id_chan_Conf2Intere_mode_read_data(avm_channel_id_chan_Conf2Intere_mode_read_data),
      // AVST avm_channel_id_chan_Conf2Intere_x_read
      .avm_channel_id_chan_Conf2Intere_x_read_valid(avm_channel_id_chan_Conf2Intere_x_read_valid),
      .avm_channel_id_chan_Conf2Intere_x_read_ready(avm_channel_id_chan_Conf2Intere_x_read_ready),
      .avm_channel_id_chan_Conf2Intere_x_read_data(avm_channel_id_chan_Conf2Intere_x_read_data),
      // AVST avm_channel_id_chan_Conf2Intere_y_read
      .avm_channel_id_chan_Conf2Intere_y_read_valid(avm_channel_id_chan_Conf2Intere_y_read_valid),
      .avm_channel_id_chan_Conf2Intere_y_read_ready(avm_channel_id_chan_Conf2Intere_y_read_ready),
      .avm_channel_id_chan_Conf2Intere_y_read_data(avm_channel_id_chan_Conf2Intere_y_read_data),
      // AVST avm_channel_id_chan_Conf2Intere_z_read
      .avm_channel_id_chan_Conf2Intere_z_read_valid(avm_channel_id_chan_Conf2Intere_z_read_valid),
      .avm_channel_id_chan_Conf2Intere_z_read_ready(avm_channel_id_chan_Conf2Intere_z_read_ready),
      .avm_channel_id_chan_Conf2Intere_z_read_data(avm_channel_id_chan_Conf2Intere_z_read_data),
      // AVST avm_channel_id_chan_Intere2Store_active_write
      .avm_channel_id_chan_Intere2Store_active_write_valid(avm_channel_id_chan_Intere2Store_active_write_valid),
      .avm_channel_id_chan_Intere2Store_active_write_ready(avm_channel_id_chan_Intere2Store_active_write_ready),
      .avm_channel_id_chan_Intere2Store_active_write_data(avm_channel_id_chan_Intere2Store_active_write_data),
      .avm_channel_id_chan_Intere2Store_active_write_almostfull(avm_channel_id_chan_Intere2Store_active_write_almostfull),
      // AVST avm_channel_id_chan_Intere2Store_cnt_write
      .avm_channel_id_chan_Intere2Store_cnt_write_valid(avm_channel_id_chan_Intere2Store_cnt_write_valid),
      .avm_channel_id_chan_Intere2Store_cnt_write_ready(avm_channel_id_chan_Intere2Store_cnt_write_ready),
      .avm_channel_id_chan_Intere2Store_cnt_write_data(avm_channel_id_chan_Intere2Store_cnt_write_data),
      .avm_channel_id_chan_Intere2Store_cnt_write_almostfull(avm_channel_id_chan_Intere2Store_cnt_write_almostfull),
      // AVST avm_channel_id_chan_Intere2Store_intere_write
      .avm_channel_id_chan_Intere2Store_intere_write_valid(avm_channel_id_chan_Intere2Store_intere_write_valid),
      .avm_channel_id_chan_Intere2Store_intere_write_ready(avm_channel_id_chan_Intere2Store_intere_write_ready),
      .avm_channel_id_chan_Intere2Store_intere_write_data(avm_channel_id_chan_Intere2Store_intere_write_data),
      .avm_channel_id_chan_Intere2Store_intere_write_almostfull(avm_channel_id_chan_Intere2Store_intere_write_almostfull),
      // AVST avm_channel_id_chan_Intere2Store_mode_write
      .avm_channel_id_chan_Intere2Store_mode_write_valid(avm_channel_id_chan_Intere2Store_mode_write_valid),
      .avm_channel_id_chan_Intere2Store_mode_write_ready(avm_channel_id_chan_Intere2Store_mode_write_ready),
      .avm_channel_id_chan_Intere2Store_mode_write_data(avm_channel_id_chan_Intere2Store_mode_write_data),
      .avm_channel_id_chan_Intere2Store_mode_write_almostfull(avm_channel_id_chan_Intere2Store_mode_write_almostfull),
      // AVM p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_InterE_enable[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_InterE_read[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_InterE_write[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_InterE_burstcount[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_InterE_address[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_InterE_writedata[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_InterE_byteenable[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_InterE_waitrequest[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_InterE_readdata[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_InterE_readdatavalid[0])
   );

   // INST Krnl_InterE_start_elem_inst_0 of acl_start_signal_chain_element
   acl_start_signal_chain_element Krnl_InterE_start_elem_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start_in(Krnl_InterE_start_chain[0]),
      .start_kernel(Krnl_InterE_start_kernel_copy[0]),
      .start_finish_detector(Krnl_InterE_start_task_fd[0]),
      .start_finish_chain_element(Krnl_InterE_start_finish_element[0]),
      .start_chain()
   );

   assign Krnl_IntraE_start_chain[0] = Krnl_IntraE_start;
   assign Krnl_IntraE_finish_chain[0] = 1'b1;
   assign Krnl_IntraE_cra_pending_write = |Krnl_IntraE_pending_write;
   assign Krnl_IntraE_cra_lsu_active = |Krnl_IntraE_lsu_active;
   assign Krnl_IntraE_cra_valid_in = |Krnl_IntraE_valid_in;
   assign Krnl_IntraE_stall_in = 0;
   // INST Krnl_IntraE_workgroup_dispatcher of acl_work_group_dispatcher
   acl_work_group_dispatcher
   #(
      .WIDTH(32),
      .NUM_COPIES(1),
      .RUN_FOREVER(0)
   )
   Krnl_IntraE_workgroup_dispatcher
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_IntraE_start),
      .num_groups(Krnl_IntraE_num_groups),
      .local_size(Krnl_IntraE_local_size),
      .stall_in(Krnl_IntraE_wg_disp_stall_in),
      .valid_out(Krnl_IntraE_wg_disp_valid_out),
      .group_id_out(Krnl_IntraE_wg_disp_group_id_out),
      .global_id_base_out(Krnl_IntraE_wg_disp_global_id_base_out),
      .start_out(Krnl_IntraE_wg_disp_start_out),
      .dispatched_all_groups(Krnl_IntraE_wg_disp_dispatched_all_groups)
   );

   // INST Krnl_IntraE_finish_detector of acl_kernel_finish_detector
   acl_kernel_finish_detector
   #(
      .NUM_COPIES(1),
      .WG_SIZE_W(32),
      .GLOBAL_ID_W(32),
      .TESSELLATION_SIZE(19)
   )
   Krnl_IntraE_finish_detector
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_IntraE_start),
      .wg_size(Krnl_IntraE_wg_size),
      .wg_dispatch_valid_out(Krnl_IntraE_wg_disp_valid_out),
      .wg_dispatch_stall_in(Krnl_IntraE_wg_disp_stall_in),
      .dispatched_all_groups(Krnl_IntraE_wg_disp_dispatched_all_groups),
      .kernel_copy_valid_out(Krnl_IntraE_valid_out),
      .kernel_copy_stall_in(Krnl_IntraE_stall_in),
      .pending_writes(Krnl_IntraE_cra_pending_write),
      .finish(Krnl_IntraE_finish)
   );

   // INST Krnl_IntraE_cra_slave_inst of Krnl_IntraE_function_cra_slave
   Krnl_IntraE_function_cra_slave Krnl_IntraE_cra_slave_inst
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_IntraE_start),
      .finish(Krnl_IntraE_finish),
      .global_offset_0(Krnl_IntraE_global_offset[0]),
      .global_offset_1(Krnl_IntraE_global_offset[1]),
      .global_offset_2(Krnl_IntraE_global_offset[2]),
      .work_dim(Krnl_IntraE_work_dim),
      .has_a_lsu_active(Krnl_IntraE_cra_lsu_active),
      .has_a_write_pending(Krnl_IntraE_cra_pending_write),
      .valid_in(Krnl_IntraE_cra_valid_in),
      .global_size_0(Krnl_IntraE_global_size[0]),
      .global_size_1(Krnl_IntraE_global_size[1]),
      .global_size_2(Krnl_IntraE_global_size[2]),
      .num_groups_0(Krnl_IntraE_num_groups[0]),
      .num_groups_1(Krnl_IntraE_num_groups[1]),
      .num_groups_2(Krnl_IntraE_num_groups[2]),
      .local_size_0(Krnl_IntraE_local_size[0]),
      .local_size_1(Krnl_IntraE_local_size[1]),
      .local_size_2(Krnl_IntraE_local_size[2]),
      .workgroup_size(Krnl_IntraE_wg_size),
      .kernel_arguments(Krnl_IntraE_kernel_arguments),
      .cra_irq(kernel_irqs[3]),
      // AVS avs_cra
      .avs_cra_enable(avs_Krnl_IntraE_cra_enable),
      .avs_cra_read(avs_Krnl_IntraE_cra_read),
      .avs_cra_write(avs_Krnl_IntraE_cra_write),
      .avs_cra_address(avs_Krnl_IntraE_cra_address),
      .avs_cra_writedata(avs_Krnl_IntraE_cra_writedata),
      .avs_cra_byteenable(avs_Krnl_IntraE_cra_byteenable),
      .avs_cra_readdata(avs_Krnl_IntraE_cra_readdata),
      .avs_cra_readdatavalid(avs_Krnl_IntraE_cra_readdatavalid),
      .acl_counter_reset(counter_reset_Krnl_IntraE),
      .acl_counter_init(counter_init_Krnl_IntraE),
      .acl_counter_limit(counter_limit_Krnl_IntraE),
      .acl_counter_size(counter_size_Krnl_IntraE),
      .acl_counter_full(counter_full_Krnl_IntraE)
   );

   // INST Krnl_IntraE_id_iter_inst_0 of acl_id_iterator
   acl_id_iterator
   #(
      .WIDTH(32),
      .LOCAL_WIDTH_X(1),
      .LOCAL_WIDTH_Y(1),
      .LOCAL_WIDTH_Z(1),
      .ENABLE_TESSELLATION(1)
   )
   Krnl_IntraE_id_iter_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_IntraE_wg_disp_start_out),
      .valid_in(Krnl_IntraE_wg_disp_valid_out[0]),
      .stall_out(Krnl_IntraE_wg_disp_stall_in[0]),
      .stall_in(Krnl_IntraE_stall_out[0]),
      .valid_out(Krnl_IntraE_valid_in[0]),
      .group_id_in(Krnl_IntraE_wg_disp_group_id_out),
      .global_id_base_in(Krnl_IntraE_wg_disp_global_id_base_out),
      .local_size(Krnl_IntraE_local_size),
      .global_size(Krnl_IntraE_global_size),
      .local_id(Krnl_IntraE_local_id[0]),
      .global_id(Krnl_IntraE_global_id[0]),
      .group_id(Krnl_IntraE_group_id[0])
   );

   // INST Krnl_IntraE_inst_0 of Krnl_IntraE_top_wrapper_0
   Krnl_IntraE_top_wrapper_0 Krnl_IntraE_inst_0
   (
      .start(Krnl_IntraE_start_kernel_copy[0]),
      .kernel_arguments(Krnl_IntraE_kernel_arguments),
      .work_dim(Krnl_IntraE_work_dim),
      .global_offset(Krnl_IntraE_global_offset),
      .kernel_valid_out(Krnl_IntraE_valid_out[0]),
      .has_a_write_pending(Krnl_IntraE_pending_write[0]),
      .has_a_lsu_active(Krnl_IntraE_lsu_active[0]),
      .global_id(Krnl_IntraE_global_id[0]),
      .local_id(Krnl_IntraE_local_id[0]),
      .group_id(Krnl_IntraE_group_id[0]),
      .global_size(Krnl_IntraE_global_size),
      .local_size(Krnl_IntraE_local_size),
      .num_groups(Krnl_IntraE_num_groups),
      .workgroup_size(Krnl_IntraE_wg_size),
      .kernel_stall_out(Krnl_IntraE_stall_out[0]),
      .kernel_valid_in(Krnl_IntraE_valid_in[0]),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_ld__inst0
      .avm_local_bb0_ld__inst0_enable(avm_kernel_rd_enable[35]),
      .avm_local_bb0_ld__inst0_read(avm_kernel_rd_read[35]),
      .avm_local_bb0_ld__inst0_write(avm_kernel_rd_write[35]),
      .avm_local_bb0_ld__inst0_burstcount(avm_kernel_rd_burstcount[35]),
      .avm_local_bb0_ld__inst0_address(avm_kernel_rd_address[35]),
      .avm_local_bb0_ld__inst0_writedata(avm_kernel_rd_writedata[35]),
      .avm_local_bb0_ld__inst0_byteenable(avm_kernel_rd_byteenable[35]),
      .avm_local_bb0_ld__inst0_waitrequest(avm_kernel_rd_waitrequest[35]),
      .avm_local_bb0_ld__inst0_readdata(avm_kernel_rd_readdata[35]),
      .avm_local_bb0_ld__inst0_readdatavalid(avm_kernel_rd_readdatavalid[35]),
      .avm_local_bb0_ld__inst0_writeack(avm_kernel_rd_writeack[35]),
      // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable(avm_kernel_rd_enable[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read(avm_kernel_rd_read[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write(avm_kernel_rd_write[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount(avm_kernel_rd_burstcount[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address(avm_kernel_rd_address[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata(avm_kernel_rd_writedata[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable(avm_kernel_rd_byteenable[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest(avm_kernel_rd_waitrequest[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata(avm_kernel_rd_readdata[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid(avm_kernel_rd_readdatavalid[36]),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack(avm_kernel_rd_writeack[36]),
      // AVM avm_local_bb3_st__inst0
      .avm_local_bb3_st__inst0_enable(avm_kernel_wr_enable[5]),
      .avm_local_bb3_st__inst0_read(avm_kernel_wr_read[5]),
      .avm_local_bb3_st__inst0_write(avm_kernel_wr_write[5]),
      .avm_local_bb3_st__inst0_burstcount(avm_kernel_wr_burstcount[5]),
      .avm_local_bb3_st__inst0_address(avm_kernel_wr_address[5]),
      .avm_local_bb3_st__inst0_writedata(avm_kernel_wr_writedata[5]),
      .avm_local_bb3_st__inst0_byteenable(avm_kernel_wr_byteenable[5]),
      .avm_local_bb3_st__inst0_waitrequest(avm_kernel_wr_waitrequest[5]),
      .avm_local_bb3_st__inst0_readdata(avm_kernel_wr_readdata[5]),
      .avm_local_bb3_st__inst0_readdatavalid(avm_kernel_wr_readdatavalid[5]),
      .avm_local_bb3_st__inst0_writeack(avm_kernel_wr_writeack[5]),
      // AVM avm_local_bb4_ld__u18_inst0
      .avm_local_bb4_ld__u18_inst0_enable(avm_kernel_rd_enable[37]),
      .avm_local_bb4_ld__u18_inst0_read(avm_kernel_rd_read[37]),
      .avm_local_bb4_ld__u18_inst0_write(avm_kernel_rd_write[37]),
      .avm_local_bb4_ld__u18_inst0_burstcount(avm_kernel_rd_burstcount[37]),
      .avm_local_bb4_ld__u18_inst0_address(avm_kernel_rd_address[37]),
      .avm_local_bb4_ld__u18_inst0_writedata(avm_kernel_rd_writedata[37]),
      .avm_local_bb4_ld__u18_inst0_byteenable(avm_kernel_rd_byteenable[37]),
      .avm_local_bb4_ld__u18_inst0_waitrequest(avm_kernel_rd_waitrequest[37]),
      .avm_local_bb4_ld__u18_inst0_readdata(avm_kernel_rd_readdata[37]),
      .avm_local_bb4_ld__u18_inst0_readdatavalid(avm_kernel_rd_readdatavalid[37]),
      .avm_local_bb4_ld__u18_inst0_writeack(avm_kernel_rd_writeack[37]),
      // AVM avm_local_bb4_ld__u19_inst0
      .avm_local_bb4_ld__u19_inst0_enable(avm_kernel_rd_enable[38]),
      .avm_local_bb4_ld__u19_inst0_read(avm_kernel_rd_read[38]),
      .avm_local_bb4_ld__u19_inst0_write(avm_kernel_rd_write[38]),
      .avm_local_bb4_ld__u19_inst0_burstcount(avm_kernel_rd_burstcount[38]),
      .avm_local_bb4_ld__u19_inst0_address(avm_kernel_rd_address[38]),
      .avm_local_bb4_ld__u19_inst0_writedata(avm_kernel_rd_writedata[38]),
      .avm_local_bb4_ld__u19_inst0_byteenable(avm_kernel_rd_byteenable[38]),
      .avm_local_bb4_ld__u19_inst0_waitrequest(avm_kernel_rd_waitrequest[38]),
      .avm_local_bb4_ld__u19_inst0_readdata(avm_kernel_rd_readdata[38]),
      .avm_local_bb4_ld__u19_inst0_readdatavalid(avm_kernel_rd_readdatavalid[38]),
      .avm_local_bb4_ld__u19_inst0_writeack(avm_kernel_rd_writeack[38]),
      // AVM avm_local_bb4_ld__u20_inst0
      .avm_local_bb4_ld__u20_inst0_enable(avm_kernel_rd_enable[39]),
      .avm_local_bb4_ld__u20_inst0_read(avm_kernel_rd_read[39]),
      .avm_local_bb4_ld__u20_inst0_write(avm_kernel_rd_write[39]),
      .avm_local_bb4_ld__u20_inst0_burstcount(avm_kernel_rd_burstcount[39]),
      .avm_local_bb4_ld__u20_inst0_address(avm_kernel_rd_address[39]),
      .avm_local_bb4_ld__u20_inst0_writedata(avm_kernel_rd_writedata[39]),
      .avm_local_bb4_ld__u20_inst0_byteenable(avm_kernel_rd_byteenable[39]),
      .avm_local_bb4_ld__u20_inst0_waitrequest(avm_kernel_rd_waitrequest[39]),
      .avm_local_bb4_ld__u20_inst0_readdata(avm_kernel_rd_readdata[39]),
      .avm_local_bb4_ld__u20_inst0_readdatavalid(avm_kernel_rd_readdatavalid[39]),
      .avm_local_bb4_ld__u20_inst0_writeack(avm_kernel_rd_writeack[39]),
      // AVM avm_local_bb4_ld__u21_inst0
      .avm_local_bb4_ld__u21_inst0_enable(avm_kernel_rd_enable[40]),
      .avm_local_bb4_ld__u21_inst0_read(avm_kernel_rd_read[40]),
      .avm_local_bb4_ld__u21_inst0_write(avm_kernel_rd_write[40]),
      .avm_local_bb4_ld__u21_inst0_burstcount(avm_kernel_rd_burstcount[40]),
      .avm_local_bb4_ld__u21_inst0_address(avm_kernel_rd_address[40]),
      .avm_local_bb4_ld__u21_inst0_writedata(avm_kernel_rd_writedata[40]),
      .avm_local_bb4_ld__u21_inst0_byteenable(avm_kernel_rd_byteenable[40]),
      .avm_local_bb4_ld__u21_inst0_waitrequest(avm_kernel_rd_waitrequest[40]),
      .avm_local_bb4_ld__u21_inst0_readdata(avm_kernel_rd_readdata[40]),
      .avm_local_bb4_ld__u21_inst0_readdatavalid(avm_kernel_rd_readdatavalid[40]),
      .avm_local_bb4_ld__u21_inst0_writeack(avm_kernel_rd_writeack[40]),
      // AVM avm_local_bb4_ld__u22_inst0
      .avm_local_bb4_ld__u22_inst0_enable(avm_kernel_rd_enable[41]),
      .avm_local_bb4_ld__u22_inst0_read(avm_kernel_rd_read[41]),
      .avm_local_bb4_ld__u22_inst0_write(avm_kernel_rd_write[41]),
      .avm_local_bb4_ld__u22_inst0_burstcount(avm_kernel_rd_burstcount[41]),
      .avm_local_bb4_ld__u22_inst0_address(avm_kernel_rd_address[41]),
      .avm_local_bb4_ld__u22_inst0_writedata(avm_kernel_rd_writedata[41]),
      .avm_local_bb4_ld__u22_inst0_byteenable(avm_kernel_rd_byteenable[41]),
      .avm_local_bb4_ld__u22_inst0_waitrequest(avm_kernel_rd_waitrequest[41]),
      .avm_local_bb4_ld__u22_inst0_readdata(avm_kernel_rd_readdata[41]),
      .avm_local_bb4_ld__u22_inst0_readdatavalid(avm_kernel_rd_readdatavalid[41]),
      .avm_local_bb4_ld__u22_inst0_writeack(avm_kernel_rd_writeack[41]),
      // AVM avm_local_bb4_ld__u23_inst0
      .avm_local_bb4_ld__u23_inst0_enable(avm_kernel_rd_enable[42]),
      .avm_local_bb4_ld__u23_inst0_read(avm_kernel_rd_read[42]),
      .avm_local_bb4_ld__u23_inst0_write(avm_kernel_rd_write[42]),
      .avm_local_bb4_ld__u23_inst0_burstcount(avm_kernel_rd_burstcount[42]),
      .avm_local_bb4_ld__u23_inst0_address(avm_kernel_rd_address[42]),
      .avm_local_bb4_ld__u23_inst0_writedata(avm_kernel_rd_writedata[42]),
      .avm_local_bb4_ld__u23_inst0_byteenable(avm_kernel_rd_byteenable[42]),
      .avm_local_bb4_ld__u23_inst0_waitrequest(avm_kernel_rd_waitrequest[42]),
      .avm_local_bb4_ld__u23_inst0_readdata(avm_kernel_rd_readdata[42]),
      .avm_local_bb4_ld__u23_inst0_readdatavalid(avm_kernel_rd_readdatavalid[42]),
      .avm_local_bb4_ld__u23_inst0_writeack(avm_kernel_rd_writeack[42]),
      // AVM avm_local_bb4_ld__u24_inst0
      .avm_local_bb4_ld__u24_inst0_enable(avm_kernel_rd_enable[43]),
      .avm_local_bb4_ld__u24_inst0_read(avm_kernel_rd_read[43]),
      .avm_local_bb4_ld__u24_inst0_write(avm_kernel_rd_write[43]),
      .avm_local_bb4_ld__u24_inst0_burstcount(avm_kernel_rd_burstcount[43]),
      .avm_local_bb4_ld__u24_inst0_address(avm_kernel_rd_address[43]),
      .avm_local_bb4_ld__u24_inst0_writedata(avm_kernel_rd_writedata[43]),
      .avm_local_bb4_ld__u24_inst0_byteenable(avm_kernel_rd_byteenable[43]),
      .avm_local_bb4_ld__u24_inst0_waitrequest(avm_kernel_rd_waitrequest[43]),
      .avm_local_bb4_ld__u24_inst0_readdata(avm_kernel_rd_readdata[43]),
      .avm_local_bb4_ld__u24_inst0_readdatavalid(avm_kernel_rd_readdatavalid[43]),
      .avm_local_bb4_ld__u24_inst0_writeack(avm_kernel_rd_writeack[43]),
      // AVM avm_local_bb4_ld__u25_inst0
      .avm_local_bb4_ld__u25_inst0_enable(avm_kernel_rd_enable[44]),
      .avm_local_bb4_ld__u25_inst0_read(avm_kernel_rd_read[44]),
      .avm_local_bb4_ld__u25_inst0_write(avm_kernel_rd_write[44]),
      .avm_local_bb4_ld__u25_inst0_burstcount(avm_kernel_rd_burstcount[44]),
      .avm_local_bb4_ld__u25_inst0_address(avm_kernel_rd_address[44]),
      .avm_local_bb4_ld__u25_inst0_writedata(avm_kernel_rd_writedata[44]),
      .avm_local_bb4_ld__u25_inst0_byteenable(avm_kernel_rd_byteenable[44]),
      .avm_local_bb4_ld__u25_inst0_waitrequest(avm_kernel_rd_waitrequest[44]),
      .avm_local_bb4_ld__u25_inst0_readdata(avm_kernel_rd_readdata[44]),
      .avm_local_bb4_ld__u25_inst0_readdatavalid(avm_kernel_rd_readdatavalid[44]),
      .avm_local_bb4_ld__u25_inst0_writeack(avm_kernel_rd_writeack[44]),
      // AVM avm_local_bb4_ld__u26_inst0
      .avm_local_bb4_ld__u26_inst0_enable(avm_kernel_rd_enable[45]),
      .avm_local_bb4_ld__u26_inst0_read(avm_kernel_rd_read[45]),
      .avm_local_bb4_ld__u26_inst0_write(avm_kernel_rd_write[45]),
      .avm_local_bb4_ld__u26_inst0_burstcount(avm_kernel_rd_burstcount[45]),
      .avm_local_bb4_ld__u26_inst0_address(avm_kernel_rd_address[45]),
      .avm_local_bb4_ld__u26_inst0_writedata(avm_kernel_rd_writedata[45]),
      .avm_local_bb4_ld__u26_inst0_byteenable(avm_kernel_rd_byteenable[45]),
      .avm_local_bb4_ld__u26_inst0_waitrequest(avm_kernel_rd_waitrequest[45]),
      .avm_local_bb4_ld__u26_inst0_readdata(avm_kernel_rd_readdata[45]),
      .avm_local_bb4_ld__u26_inst0_readdatavalid(avm_kernel_rd_readdatavalid[45]),
      .avm_local_bb4_ld__u26_inst0_writeack(avm_kernel_rd_writeack[45]),
      // AVM avm_local_bb4_ld__u27_inst0
      .avm_local_bb4_ld__u27_inst0_enable(avm_kernel_rd_enable[46]),
      .avm_local_bb4_ld__u27_inst0_read(avm_kernel_rd_read[46]),
      .avm_local_bb4_ld__u27_inst0_write(avm_kernel_rd_write[46]),
      .avm_local_bb4_ld__u27_inst0_burstcount(avm_kernel_rd_burstcount[46]),
      .avm_local_bb4_ld__u27_inst0_address(avm_kernel_rd_address[46]),
      .avm_local_bb4_ld__u27_inst0_writedata(avm_kernel_rd_writedata[46]),
      .avm_local_bb4_ld__u27_inst0_byteenable(avm_kernel_rd_byteenable[46]),
      .avm_local_bb4_ld__u27_inst0_waitrequest(avm_kernel_rd_waitrequest[46]),
      .avm_local_bb4_ld__u27_inst0_readdata(avm_kernel_rd_readdata[46]),
      .avm_local_bb4_ld__u27_inst0_readdatavalid(avm_kernel_rd_readdatavalid[46]),
      .avm_local_bb4_ld__u27_inst0_writeack(avm_kernel_rd_writeack[46]),
      // AVM avm_local_bb4_ld__u28_inst0
      .avm_local_bb4_ld__u28_inst0_enable(avm_kernel_rd_enable[47]),
      .avm_local_bb4_ld__u28_inst0_read(avm_kernel_rd_read[47]),
      .avm_local_bb4_ld__u28_inst0_write(avm_kernel_rd_write[47]),
      .avm_local_bb4_ld__u28_inst0_burstcount(avm_kernel_rd_burstcount[47]),
      .avm_local_bb4_ld__u28_inst0_address(avm_kernel_rd_address[47]),
      .avm_local_bb4_ld__u28_inst0_writedata(avm_kernel_rd_writedata[47]),
      .avm_local_bb4_ld__u28_inst0_byteenable(avm_kernel_rd_byteenable[47]),
      .avm_local_bb4_ld__u28_inst0_waitrequest(avm_kernel_rd_waitrequest[47]),
      .avm_local_bb4_ld__u28_inst0_readdata(avm_kernel_rd_readdata[47]),
      .avm_local_bb4_ld__u28_inst0_readdatavalid(avm_kernel_rd_readdatavalid[47]),
      .avm_local_bb4_ld__u28_inst0_writeack(avm_kernel_rd_writeack[47]),
      // AVM avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_enable(avm_kernel_rd_enable[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_read(avm_kernel_rd_read[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_write(avm_kernel_rd_write[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_burstcount(avm_kernel_rd_burstcount[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_address(avm_kernel_rd_address[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writedata(avm_kernel_rd_writedata[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_byteenable(avm_kernel_rd_byteenable[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_waitrequest(avm_kernel_rd_waitrequest[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdata(avm_kernel_rd_readdata[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid(avm_kernel_rd_readdatavalid[48]),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writeack(avm_kernel_rd_writeack[48]),
      // AVST avm_channel_id_chan_Conf2Intrae_active_read
      .avm_channel_id_chan_Conf2Intrae_active_read_valid(avm_channel_id_chan_Conf2Intrae_active_read_valid),
      .avm_channel_id_chan_Conf2Intrae_active_read_ready(avm_channel_id_chan_Conf2Intrae_active_read_ready),
      .avm_channel_id_chan_Conf2Intrae_active_read_data(avm_channel_id_chan_Conf2Intrae_active_read_data),
      // AVST avm_channel_id_chan_Conf2Intrae_cnt_read
      .avm_channel_id_chan_Conf2Intrae_cnt_read_valid(avm_channel_id_chan_Conf2Intrae_cnt_read_valid),
      .avm_channel_id_chan_Conf2Intrae_cnt_read_ready(avm_channel_id_chan_Conf2Intrae_cnt_read_ready),
      .avm_channel_id_chan_Conf2Intrae_cnt_read_data(avm_channel_id_chan_Conf2Intrae_cnt_read_data),
      // AVST avm_channel_id_chan_Conf2Intrae_mode_read
      .avm_channel_id_chan_Conf2Intrae_mode_read_valid(avm_channel_id_chan_Conf2Intrae_mode_read_valid),
      .avm_channel_id_chan_Conf2Intrae_mode_read_ready(avm_channel_id_chan_Conf2Intrae_mode_read_ready),
      .avm_channel_id_chan_Conf2Intrae_mode_read_data(avm_channel_id_chan_Conf2Intrae_mode_read_data),
      // AVST avm_channel_id_chan_Conf2Intrae_x_read
      .avm_channel_id_chan_Conf2Intrae_x_read_valid(avm_channel_id_chan_Conf2Intrae_x_read_valid),
      .avm_channel_id_chan_Conf2Intrae_x_read_ready(avm_channel_id_chan_Conf2Intrae_x_read_ready),
      .avm_channel_id_chan_Conf2Intrae_x_read_data(avm_channel_id_chan_Conf2Intrae_x_read_data),
      // AVST avm_channel_id_chan_Conf2Intrae_y_read
      .avm_channel_id_chan_Conf2Intrae_y_read_valid(avm_channel_id_chan_Conf2Intrae_y_read_valid),
      .avm_channel_id_chan_Conf2Intrae_y_read_ready(avm_channel_id_chan_Conf2Intrae_y_read_ready),
      .avm_channel_id_chan_Conf2Intrae_y_read_data(avm_channel_id_chan_Conf2Intrae_y_read_data),
      // AVST avm_channel_id_chan_Conf2Intrae_z_read
      .avm_channel_id_chan_Conf2Intrae_z_read_valid(avm_channel_id_chan_Conf2Intrae_z_read_valid),
      .avm_channel_id_chan_Conf2Intrae_z_read_ready(avm_channel_id_chan_Conf2Intrae_z_read_ready),
      .avm_channel_id_chan_Conf2Intrae_z_read_data(avm_channel_id_chan_Conf2Intrae_z_read_data),
      // AVST avm_channel_id_chan_Intrae2Store_active_write
      .avm_channel_id_chan_Intrae2Store_active_write_valid(avm_channel_id_chan_Intrae2Store_active_write_valid),
      .avm_channel_id_chan_Intrae2Store_active_write_ready(avm_channel_id_chan_Intrae2Store_active_write_ready),
      .avm_channel_id_chan_Intrae2Store_active_write_data(avm_channel_id_chan_Intrae2Store_active_write_data),
      .avm_channel_id_chan_Intrae2Store_active_write_almostfull(avm_channel_id_chan_Intrae2Store_active_write_almostfull),
      // AVST avm_channel_id_chan_Intrae2Store_cnt_write
      .avm_channel_id_chan_Intrae2Store_cnt_write_valid(avm_channel_id_chan_Intrae2Store_cnt_write_valid),
      .avm_channel_id_chan_Intrae2Store_cnt_write_ready(avm_channel_id_chan_Intrae2Store_cnt_write_ready),
      .avm_channel_id_chan_Intrae2Store_cnt_write_data(avm_channel_id_chan_Intrae2Store_cnt_write_data),
      .avm_channel_id_chan_Intrae2Store_cnt_write_almostfull(avm_channel_id_chan_Intrae2Store_cnt_write_almostfull),
      // AVST avm_channel_id_chan_Intrae2Store_intrae_write
      .avm_channel_id_chan_Intrae2Store_intrae_write_valid(avm_channel_id_chan_Intrae2Store_intrae_write_valid),
      .avm_channel_id_chan_Intrae2Store_intrae_write_ready(avm_channel_id_chan_Intrae2Store_intrae_write_ready),
      .avm_channel_id_chan_Intrae2Store_intrae_write_data(avm_channel_id_chan_Intrae2Store_intrae_write_data),
      .avm_channel_id_chan_Intrae2Store_intrae_write_almostfull(avm_channel_id_chan_Intrae2Store_intrae_write_almostfull),
      // AVST avm_channel_id_chan_Intrae2Store_mode_write
      .avm_channel_id_chan_Intrae2Store_mode_write_valid(avm_channel_id_chan_Intrae2Store_mode_write_valid),
      .avm_channel_id_chan_Intrae2Store_mode_write_ready(avm_channel_id_chan_Intrae2Store_mode_write_ready),
      .avm_channel_id_chan_Intrae2Store_mode_write_data(avm_channel_id_chan_Intrae2Store_mode_write_data),
      .avm_channel_id_chan_Intrae2Store_mode_write_almostfull(avm_channel_id_chan_Intrae2Store_mode_write_almostfull),
      // AVM p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_IntraE_enable[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_IntraE_read[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_IntraE_write[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_IntraE_burstcount[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_IntraE_address[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_IntraE_writedata[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_IntraE_byteenable[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_IntraE_waitrequest[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_IntraE_readdata[0]),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_IntraE_readdatavalid[0])
   );

   // INST Krnl_IntraE_start_elem_inst_0 of acl_start_signal_chain_element
   acl_start_signal_chain_element Krnl_IntraE_start_elem_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start_in(Krnl_IntraE_start_chain[0]),
      .start_kernel(Krnl_IntraE_start_kernel_copy[0]),
      .start_finish_detector(Krnl_IntraE_start_task_fd[0]),
      .start_finish_chain_element(Krnl_IntraE_start_finish_element[0]),
      .start_chain()
   );

   assign Krnl_Store_start_chain[0] = Krnl_Store_start;
   assign Krnl_Store_finish_chain[0] = 1'b1;
   assign Krnl_Store_cra_pending_write = |Krnl_Store_pending_write;
   assign Krnl_Store_cra_lsu_active = |Krnl_Store_lsu_active;
   assign Krnl_Store_cra_valid_in = |Krnl_Store_valid_in;
   assign Krnl_Store_stall_in = 0;
   // INST Krnl_Store_workgroup_dispatcher of acl_work_group_dispatcher
   acl_work_group_dispatcher
   #(
      .WIDTH(32),
      .NUM_COPIES(1),
      .RUN_FOREVER(0)
   )
   Krnl_Store_workgroup_dispatcher
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Store_start),
      .num_groups(Krnl_Store_num_groups),
      .local_size(Krnl_Store_local_size),
      .stall_in(Krnl_Store_wg_disp_stall_in),
      .valid_out(Krnl_Store_wg_disp_valid_out),
      .group_id_out(Krnl_Store_wg_disp_group_id_out),
      .global_id_base_out(Krnl_Store_wg_disp_global_id_base_out),
      .start_out(Krnl_Store_wg_disp_start_out),
      .dispatched_all_groups(Krnl_Store_wg_disp_dispatched_all_groups)
   );

   // INST Krnl_Store_finish_detector of acl_kernel_finish_detector
   acl_kernel_finish_detector
   #(
      .NUM_COPIES(1),
      .WG_SIZE_W(32),
      .GLOBAL_ID_W(32),
      .TESSELLATION_SIZE(19)
   )
   Krnl_Store_finish_detector
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Store_start),
      .wg_size(Krnl_Store_wg_size),
      .wg_dispatch_valid_out(Krnl_Store_wg_disp_valid_out),
      .wg_dispatch_stall_in(Krnl_Store_wg_disp_stall_in),
      .dispatched_all_groups(Krnl_Store_wg_disp_dispatched_all_groups),
      .kernel_copy_valid_out(Krnl_Store_valid_out),
      .kernel_copy_stall_in(Krnl_Store_stall_in),
      .pending_writes(Krnl_Store_cra_pending_write),
      .finish(Krnl_Store_finish)
   );

   // INST Krnl_Store_cra_slave_inst of Krnl_Store_function_cra_slave
   Krnl_Store_function_cra_slave Krnl_Store_cra_slave_inst
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Store_start),
      .finish(Krnl_Store_finish),
      .global_offset_0(Krnl_Store_global_offset[0]),
      .global_offset_1(Krnl_Store_global_offset[1]),
      .global_offset_2(Krnl_Store_global_offset[2]),
      .work_dim(Krnl_Store_work_dim),
      .has_a_lsu_active(Krnl_Store_cra_lsu_active),
      .has_a_write_pending(Krnl_Store_cra_pending_write),
      .valid_in(Krnl_Store_cra_valid_in),
      .global_size_0(Krnl_Store_global_size[0]),
      .global_size_1(Krnl_Store_global_size[1]),
      .global_size_2(Krnl_Store_global_size[2]),
      .num_groups_0(Krnl_Store_num_groups[0]),
      .num_groups_1(Krnl_Store_num_groups[1]),
      .num_groups_2(Krnl_Store_num_groups[2]),
      .local_size_0(Krnl_Store_local_size[0]),
      .local_size_1(Krnl_Store_local_size[1]),
      .local_size_2(Krnl_Store_local_size[2]),
      .workgroup_size(Krnl_Store_wg_size),
      .kernel_arguments(Krnl_Store_kernel_arguments),
      .cra_irq(kernel_irqs[4]),
      // AVS avs_cra
      .avs_cra_enable(avs_Krnl_Store_cra_enable),
      .avs_cra_read(avs_Krnl_Store_cra_read),
      .avs_cra_write(avs_Krnl_Store_cra_write),
      .avs_cra_address(avs_Krnl_Store_cra_address),
      .avs_cra_writedata(avs_Krnl_Store_cra_writedata),
      .avs_cra_byteenable(avs_Krnl_Store_cra_byteenable),
      .avs_cra_readdata(avs_Krnl_Store_cra_readdata),
      .avs_cra_readdatavalid(avs_Krnl_Store_cra_readdatavalid),
      .acl_counter_reset(counter_reset_Krnl_Store),
      .acl_counter_init(counter_init_Krnl_Store),
      .acl_counter_limit(counter_limit_Krnl_Store),
      .acl_counter_size(counter_size_Krnl_Store),
      .acl_counter_full(counter_full_Krnl_Store)
   );

   // INST Krnl_Store_id_iter_inst_0 of acl_id_iterator
   acl_id_iterator
   #(
      .WIDTH(32),
      .LOCAL_WIDTH_X(32),
      .LOCAL_WIDTH_Y(32),
      .LOCAL_WIDTH_Z(32),
      .ENABLE_TESSELLATION(1)
   )
   Krnl_Store_id_iter_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start(Krnl_Store_wg_disp_start_out),
      .valid_in(Krnl_Store_wg_disp_valid_out[0]),
      .stall_out(Krnl_Store_wg_disp_stall_in[0]),
      .stall_in(Krnl_Store_stall_out[0]),
      .valid_out(Krnl_Store_valid_in[0]),
      .group_id_in(Krnl_Store_wg_disp_group_id_out),
      .global_id_base_in(Krnl_Store_wg_disp_global_id_base_out),
      .local_size(Krnl_Store_local_size),
      .global_size(Krnl_Store_global_size),
      .local_id(Krnl_Store_local_id[0]),
      .global_id(Krnl_Store_global_id[0]),
      .group_id(Krnl_Store_group_id[0])
   );

   // INST Krnl_Store_inst_0 of Krnl_Store_top_wrapper_0
   Krnl_Store_top_wrapper_0 Krnl_Store_inst_0
   (
      .start(Krnl_Store_start_kernel_copy[0]),
      .kernel_arguments(Krnl_Store_kernel_arguments),
      .work_dim(Krnl_Store_work_dim),
      .global_offset(Krnl_Store_global_offset),
      .kernel_valid_out(Krnl_Store_valid_out[0]),
      .has_a_write_pending(Krnl_Store_pending_write[0]),
      .has_a_lsu_active(Krnl_Store_lsu_active[0]),
      .global_id(Krnl_Store_global_id[0]),
      .local_id(Krnl_Store_local_id[0]),
      .group_id(Krnl_Store_group_id[0]),
      .global_size(Krnl_Store_global_size),
      .local_size(Krnl_Store_local_size),
      .num_groups(Krnl_Store_num_groups),
      .workgroup_size(Krnl_Store_wg_size),
      .kernel_stall_out(Krnl_Store_stall_out[0]),
      .kernel_valid_in(Krnl_Store_valid_in[0]),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb1_st__inst0
      .avm_local_bb1_st__inst0_enable(avm_kernel_wr_enable[6]),
      .avm_local_bb1_st__inst0_read(avm_kernel_wr_read[6]),
      .avm_local_bb1_st__inst0_write(avm_kernel_wr_write[6]),
      .avm_local_bb1_st__inst0_burstcount(avm_kernel_wr_burstcount[6]),
      .avm_local_bb1_st__inst0_address(avm_kernel_wr_address[6]),
      .avm_local_bb1_st__inst0_writedata(avm_kernel_wr_writedata[6]),
      .avm_local_bb1_st__inst0_byteenable(avm_kernel_wr_byteenable[6]),
      .avm_local_bb1_st__inst0_waitrequest(avm_kernel_wr_waitrequest[6]),
      .avm_local_bb1_st__inst0_readdata(avm_kernel_wr_readdata[6]),
      .avm_local_bb1_st__inst0_readdatavalid(avm_kernel_wr_readdatavalid[6]),
      .avm_local_bb1_st__inst0_writeack(avm_kernel_wr_writeack[6]),
      // AVM avm_local_bb1_st__u5_inst0
      .avm_local_bb1_st__u5_inst0_enable(avm_kernel_wr_enable[7]),
      .avm_local_bb1_st__u5_inst0_read(avm_kernel_wr_read[7]),
      .avm_local_bb1_st__u5_inst0_write(avm_kernel_wr_write[7]),
      .avm_local_bb1_st__u5_inst0_burstcount(avm_kernel_wr_burstcount[7]),
      .avm_local_bb1_st__u5_inst0_address(avm_kernel_wr_address[7]),
      .avm_local_bb1_st__u5_inst0_writedata(avm_kernel_wr_writedata[7]),
      .avm_local_bb1_st__u5_inst0_byteenable(avm_kernel_wr_byteenable[7]),
      .avm_local_bb1_st__u5_inst0_waitrequest(avm_kernel_wr_waitrequest[7]),
      .avm_local_bb1_st__u5_inst0_readdata(avm_kernel_wr_readdata[7]),
      .avm_local_bb1_st__u5_inst0_readdatavalid(avm_kernel_wr_readdatavalid[7]),
      .avm_local_bb1_st__u5_inst0_writeack(avm_kernel_wr_writeack[7]),
      // AVM avm_local_bb1_st__u7_inst0
      .avm_local_bb1_st__u7_inst0_enable(avm_kernel_wr_enable[8]),
      .avm_local_bb1_st__u7_inst0_read(avm_kernel_wr_read[8]),
      .avm_local_bb1_st__u7_inst0_write(avm_kernel_wr_write[8]),
      .avm_local_bb1_st__u7_inst0_burstcount(avm_kernel_wr_burstcount[8]),
      .avm_local_bb1_st__u7_inst0_address(avm_kernel_wr_address[8]),
      .avm_local_bb1_st__u7_inst0_writedata(avm_kernel_wr_writedata[8]),
      .avm_local_bb1_st__u7_inst0_byteenable(avm_kernel_wr_byteenable[8]),
      .avm_local_bb1_st__u7_inst0_waitrequest(avm_kernel_wr_waitrequest[8]),
      .avm_local_bb1_st__u7_inst0_readdata(avm_kernel_wr_readdata[8]),
      .avm_local_bb1_st__u7_inst0_readdatavalid(avm_kernel_wr_readdatavalid[8]),
      .avm_local_bb1_st__u7_inst0_writeack(avm_kernel_wr_writeack[8]),
      // AVM avm_local_bb1_st__u9_inst0
      .avm_local_bb1_st__u9_inst0_enable(avm_kernel_wr_enable[9]),
      .avm_local_bb1_st__u9_inst0_read(avm_kernel_wr_read[9]),
      .avm_local_bb1_st__u9_inst0_write(avm_kernel_wr_write[9]),
      .avm_local_bb1_st__u9_inst0_burstcount(avm_kernel_wr_burstcount[9]),
      .avm_local_bb1_st__u9_inst0_address(avm_kernel_wr_address[9]),
      .avm_local_bb1_st__u9_inst0_writedata(avm_kernel_wr_writedata[9]),
      .avm_local_bb1_st__u9_inst0_byteenable(avm_kernel_wr_byteenable[9]),
      .avm_local_bb1_st__u9_inst0_waitrequest(avm_kernel_wr_waitrequest[9]),
      .avm_local_bb1_st__u9_inst0_readdata(avm_kernel_wr_readdata[9]),
      .avm_local_bb1_st__u9_inst0_readdatavalid(avm_kernel_wr_readdatavalid[9]),
      .avm_local_bb1_st__u9_inst0_writeack(avm_kernel_wr_writeack[9]),
      // AVST avm_channel_id_chan_Intere2Store_active_read
      .avm_channel_id_chan_Intere2Store_active_read_valid(avm_channel_id_chan_Intere2Store_active_read_valid),
      .avm_channel_id_chan_Intere2Store_active_read_ready(avm_channel_id_chan_Intere2Store_active_read_ready),
      .avm_channel_id_chan_Intere2Store_active_read_data(avm_channel_id_chan_Intere2Store_active_read_data),
      // AVST avm_channel_id_chan_Intere2Store_cnt_read
      .avm_channel_id_chan_Intere2Store_cnt_read_valid(avm_channel_id_chan_Intere2Store_cnt_read_valid),
      .avm_channel_id_chan_Intere2Store_cnt_read_ready(avm_channel_id_chan_Intere2Store_cnt_read_ready),
      .avm_channel_id_chan_Intere2Store_cnt_read_data(avm_channel_id_chan_Intere2Store_cnt_read_data),
      // AVST avm_channel_id_chan_Intere2Store_intere_read
      .avm_channel_id_chan_Intere2Store_intere_read_valid(avm_channel_id_chan_Intere2Store_intere_read_valid),
      .avm_channel_id_chan_Intere2Store_intere_read_ready(avm_channel_id_chan_Intere2Store_intere_read_ready),
      .avm_channel_id_chan_Intere2Store_intere_read_data(avm_channel_id_chan_Intere2Store_intere_read_data),
      // AVST avm_channel_id_chan_Intere2Store_mode_read
      .avm_channel_id_chan_Intere2Store_mode_read_valid(avm_channel_id_chan_Intere2Store_mode_read_valid),
      .avm_channel_id_chan_Intere2Store_mode_read_ready(avm_channel_id_chan_Intere2Store_mode_read_ready),
      .avm_channel_id_chan_Intere2Store_mode_read_data(avm_channel_id_chan_Intere2Store_mode_read_data),
      // AVST avm_channel_id_chan_Intrae2Store_active_read
      .avm_channel_id_chan_Intrae2Store_active_read_valid(avm_channel_id_chan_Intrae2Store_active_read_valid),
      .avm_channel_id_chan_Intrae2Store_active_read_ready(avm_channel_id_chan_Intrae2Store_active_read_ready),
      .avm_channel_id_chan_Intrae2Store_active_read_data(avm_channel_id_chan_Intrae2Store_active_read_data),
      // AVST avm_channel_id_chan_Intrae2Store_cnt_read
      .avm_channel_id_chan_Intrae2Store_cnt_read_valid(avm_channel_id_chan_Intrae2Store_cnt_read_valid),
      .avm_channel_id_chan_Intrae2Store_cnt_read_ready(avm_channel_id_chan_Intrae2Store_cnt_read_ready),
      .avm_channel_id_chan_Intrae2Store_cnt_read_data(avm_channel_id_chan_Intrae2Store_cnt_read_data),
      // AVST avm_channel_id_chan_Intrae2Store_intrae_read
      .avm_channel_id_chan_Intrae2Store_intrae_read_valid(avm_channel_id_chan_Intrae2Store_intrae_read_valid),
      .avm_channel_id_chan_Intrae2Store_intrae_read_ready(avm_channel_id_chan_Intrae2Store_intrae_read_ready),
      .avm_channel_id_chan_Intrae2Store_intrae_read_data(avm_channel_id_chan_Intrae2Store_intrae_read_data),
      // AVST avm_channel_id_chan_Intrae2Store_mode_read
      .avm_channel_id_chan_Intrae2Store_mode_read_valid(avm_channel_id_chan_Intrae2Store_mode_read_valid),
      .avm_channel_id_chan_Intrae2Store_mode_read_ready(avm_channel_id_chan_Intrae2Store_mode_read_ready),
      .avm_channel_id_chan_Intrae2Store_mode_read_data(avm_channel_id_chan_Intrae2Store_mode_read_data),
      // AVM p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_Store_enable[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_Store_read[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_Store_write[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_Store_burstcount[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_Store_address[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_Store_writedata[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_Store_byteenable[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_Store_waitrequest[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_Store_readdata[0]),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_Store_readdatavalid[0]),
      // AVM p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_Store_enable[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_Store_read[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_Store_write[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_Store_burstcount[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_Store_address[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_Store_writedata[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_Store_byteenable[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_Store_waitrequest[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_Store_readdata[1]),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_Store_readdatavalid[1]),
      // AVM p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_Store_enable[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_Store_read[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_Store_write[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_Store_burstcount[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_Store_address[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_Store_writedata[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_Store_byteenable[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_Store_waitrequest[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_Store_readdata[2]),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_Store_readdatavalid[2]),
      // AVM p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_enable(printf_avm_Krnl_Store_enable[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_read(printf_avm_Krnl_Store_read[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_write(printf_avm_Krnl_Store_write[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_burstcount(printf_avm_Krnl_Store_burstcount[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_address(printf_avm_Krnl_Store_address[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_writedata(printf_avm_Krnl_Store_writedata[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_byteenable(printf_avm_Krnl_Store_byteenable[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_waitrequest(printf_avm_Krnl_Store_waitrequest[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdata(printf_avm_Krnl_Store_readdata[3]),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdatavalid(printf_avm_Krnl_Store_readdatavalid[3])
   );

   // INST Krnl_Store_start_elem_inst_0 of acl_start_signal_chain_element
   acl_start_signal_chain_element Krnl_Store_start_elem_inst_0
   (
      .clock(clock),
      .resetn(resetn),
      .start_in(Krnl_Store_start_chain[0]),
      .start_kernel(Krnl_Store_start_kernel_copy[0]),
      .start_finish_detector(Krnl_Store_start_task_fd[0]),
      .start_finish_chain_element(Krnl_Store_start_finish_element[0]),
      .start_chain()
   );

   assign kernel_irq = |kernel_irqs;
   // INST lsu_ic_top of lsu_ic_top
   lsu_ic_top
   #(
      .AWIDTH(31),
      .SHIFT(31),
      .MWIDTH_BYTES(64),
      .BURST_CNT_W(5),
      .NUM_RD_PORT(49),
      .NUM_WR_PORT(10),
      .NUM_DIMM(1),
      .ENABLE_DUAL_RING(1),
      .ENABLE_MULTIPLE_WR_RING(0),
      .ENABLE_LAST_WAIT(0),
      .ENABLE_REORDER(0),
      .NUM_REORDER(1)
   )
   lsu_ic_top
   (
      .clk(clock),
      .resetn(resetn),
      .i_rd_request(avm_kernel_rd_read),
      .i_rd_address(avm_kernel_rd_address),
      .i_rd_burstcount(avm_kernel_rd_burstcount),
      .i_wr_byteenable(avm_kernel_wr_byteenable),
      .i_wr_address(avm_kernel_wr_address),
      .i_wr_request(avm_kernel_wr_write),
      .i_wr_burstcount(avm_kernel_wr_burstcount),
      .i_wr_writedata(avm_kernel_wr_writedata),
      .i_avm_waitrequest(ic_avm_waitrequest),
      .i_avm_readdata(ic_avm_readdata),
      .i_avm_readdatavalid(ic_avm_readdatavalid),
      .o_avm_byteenable(ic_avm_byteenable),
      .o_avm_address(ic_avm_address),
      .o_avm_read(ic_avm_read),
      .o_avm_write(ic_avm_write),
      .o_avm_burstcount(ic_avm_burstcount),
      .o_wr_waitrequest(avm_kernel_wr_waitrequest),
      .o_avm_writedata(ic_avm_writedata),
      .o_avm_writeack(avm_kernel_wr_writeack),
      .o_rd_waitrequest(avm_kernel_rd_waitrequest),
      .o_avm_readdata(avm_kernel_rd_readdata),
      .o_avm_readdatavalid(avm_kernel_rd_readdatavalid)
   );

   assign avm_memgmem0_DDR_port_0_0_rw_read = ic_avm_read[0];
   assign avm_memgmem0_DDR_port_0_0_rw_write = ic_avm_write[0];
   assign avm_memgmem0_DDR_port_0_0_rw_burstcount = ic_avm_burstcount[0];
   assign avm_memgmem0_DDR_port_0_0_rw_address = ic_avm_address[0];
   assign avm_memgmem0_DDR_port_0_0_rw_writedata = ic_avm_writedata[0];
   assign avm_memgmem0_DDR_port_0_0_rw_byteenable = ic_avm_byteenable[0];
   assign ic_avm_waitrequest[0] = avm_memgmem0_DDR_port_0_0_rw_waitrequest;
   assign ic_avm_readdata[0] = avm_memgmem0_DDR_port_0_0_rw_readdata;
   assign ic_avm_readdatavalid[0] = avm_memgmem0_DDR_port_0_0_rw_readdatavalid;
   assign avs_channel_id_chan_GA2Conf_genotype_write_valid = avm_channel_id_chan_GA2Conf_genotype_write_valid;
   assign avm_channel_id_chan_GA2Conf_genotype_write_ready = avs_channel_id_chan_GA2Conf_genotype_write_ready;
   assign avs_channel_id_chan_GA2Conf_genotype_write_data = avm_channel_id_chan_GA2Conf_genotype_write_data;
   assign avm_channel_id_chan_GA2Conf_genotype_read_valid = avs_channel_id_chan_GA2Conf_genotype_read_valid;
   assign avs_channel_id_chan_GA2Conf_genotype_read_ready = avm_channel_id_chan_GA2Conf_genotype_read_ready;
   assign avm_channel_id_chan_GA2Conf_genotype_read_data = avs_channel_id_chan_GA2Conf_genotype_read_data;
   // INST channel_id_chan_GA2Conf_genotype_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(39),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_GA2Conf_genotype_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_GA2Conf_genotype_write_valid),
      .avst_in_ready(avs_channel_id_chan_GA2Conf_genotype_write_ready),
      .avst_in_data(avs_channel_id_chan_GA2Conf_genotype_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_GA2Conf_genotype_read_valid),
      .avst_out_ready(avs_channel_id_chan_GA2Conf_genotype_read_ready),
      .avst_out_data(avs_channel_id_chan_GA2Conf_genotype_read_data),
      .profile_fifosize(avs_channel_id_chan_GA2Conf_genotype_fifosize),
      .almost_full(avs_channel_id_chan_GA2Conf_genotype_write_almostfull)
   );

   assign avm_channel_id_chan_GA2Conf_genotype_write_almostfull = avs_channel_id_chan_GA2Conf_genotype_write_almostfull;
   assign avs_channel_id_chan_GA2Conf_active_write_valid = avm_channel_id_chan_GA2Conf_active_write_valid;
   assign avm_channel_id_chan_GA2Conf_active_write_ready = avs_channel_id_chan_GA2Conf_active_write_ready;
   assign avs_channel_id_chan_GA2Conf_active_write_data = avm_channel_id_chan_GA2Conf_active_write_data;
   assign avm_channel_id_chan_GA2Conf_active_read_valid = avs_channel_id_chan_GA2Conf_active_read_valid;
   assign avs_channel_id_chan_GA2Conf_active_read_ready = avm_channel_id_chan_GA2Conf_active_read_ready;
   assign avm_channel_id_chan_GA2Conf_active_read_data = avs_channel_id_chan_GA2Conf_active_read_data;
   // INST channel_id_chan_GA2Conf_active_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(48),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(0),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_GA2Conf_active_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_GA2Conf_active_write_valid),
      .avst_in_ready(avs_channel_id_chan_GA2Conf_active_write_ready),
      .avst_in_data(avs_channel_id_chan_GA2Conf_active_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_GA2Conf_active_read_valid),
      .avst_out_ready(avs_channel_id_chan_GA2Conf_active_read_ready),
      .avst_out_data(avs_channel_id_chan_GA2Conf_active_read_data),
      .profile_fifosize(avs_channel_id_chan_GA2Conf_active_fifosize),
      .almost_full(avs_channel_id_chan_GA2Conf_active_write_almostfull)
   );

   assign avm_channel_id_chan_GA2Conf_active_write_almostfull = avs_channel_id_chan_GA2Conf_active_write_almostfull;
   assign avs_channel_id_chan_GA2Conf_cnt_write_valid = avm_channel_id_chan_GA2Conf_cnt_write_valid;
   assign avm_channel_id_chan_GA2Conf_cnt_write_ready = avs_channel_id_chan_GA2Conf_cnt_write_ready;
   assign avs_channel_id_chan_GA2Conf_cnt_write_data = avm_channel_id_chan_GA2Conf_cnt_write_data;
   assign avm_channel_id_chan_GA2Conf_cnt_read_valid = avs_channel_id_chan_GA2Conf_cnt_read_valid;
   assign avs_channel_id_chan_GA2Conf_cnt_read_ready = avm_channel_id_chan_GA2Conf_cnt_read_ready;
   assign avm_channel_id_chan_GA2Conf_cnt_read_data = avs_channel_id_chan_GA2Conf_cnt_read_data;
   // INST channel_id_chan_GA2Conf_cnt_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(0),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(0),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_GA2Conf_cnt_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_GA2Conf_cnt_write_valid),
      .avst_in_ready(avs_channel_id_chan_GA2Conf_cnt_write_ready),
      .avst_in_data(avs_channel_id_chan_GA2Conf_cnt_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_GA2Conf_cnt_read_valid),
      .avst_out_ready(avs_channel_id_chan_GA2Conf_cnt_read_ready),
      .avst_out_data(avs_channel_id_chan_GA2Conf_cnt_read_data),
      .profile_fifosize(avs_channel_id_chan_GA2Conf_cnt_fifosize),
      .almost_full(avs_channel_id_chan_GA2Conf_cnt_write_almostfull)
   );

   assign avm_channel_id_chan_GA2Conf_cnt_write_almostfull = avs_channel_id_chan_GA2Conf_cnt_write_almostfull;
   assign avs_channel_id_chan_GA2Conf_mode_write_valid = avm_channel_id_chan_GA2Conf_mode_write_valid;
   assign avm_channel_id_chan_GA2Conf_mode_write_ready = avs_channel_id_chan_GA2Conf_mode_write_ready;
   assign avs_channel_id_chan_GA2Conf_mode_write_data = avm_channel_id_chan_GA2Conf_mode_write_data;
   assign avm_channel_id_chan_GA2Conf_mode_read_valid = avs_channel_id_chan_GA2Conf_mode_read_valid;
   assign avs_channel_id_chan_GA2Conf_mode_read_ready = avm_channel_id_chan_GA2Conf_mode_read_ready;
   assign avm_channel_id_chan_GA2Conf_mode_read_data = avs_channel_id_chan_GA2Conf_mode_read_data;
   // INST channel_id_chan_GA2Conf_mode_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(0),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(0),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_GA2Conf_mode_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_GA2Conf_mode_write_valid),
      .avst_in_ready(avs_channel_id_chan_GA2Conf_mode_write_ready),
      .avst_in_data(avs_channel_id_chan_GA2Conf_mode_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_GA2Conf_mode_read_valid),
      .avst_out_ready(avs_channel_id_chan_GA2Conf_mode_read_ready),
      .avst_out_data(avs_channel_id_chan_GA2Conf_mode_read_data),
      .profile_fifosize(avs_channel_id_chan_GA2Conf_mode_fifosize),
      .almost_full(avs_channel_id_chan_GA2Conf_mode_write_almostfull)
   );

   assign avm_channel_id_chan_GA2Conf_mode_write_almostfull = avs_channel_id_chan_GA2Conf_mode_write_almostfull;
   assign avs_channel_id_chan_Conf2Intere_active_write_valid = avm_channel_id_chan_Conf2Intere_active_write_valid;
   assign avm_channel_id_chan_Conf2Intere_active_write_ready = avs_channel_id_chan_Conf2Intere_active_write_ready;
   assign avs_channel_id_chan_Conf2Intere_active_write_data = avm_channel_id_chan_Conf2Intere_active_write_data;
   assign avm_channel_id_chan_Conf2Intere_active_read_valid = avs_channel_id_chan_Conf2Intere_active_read_valid;
   assign avs_channel_id_chan_Conf2Intere_active_read_ready = avm_channel_id_chan_Conf2Intere_active_read_ready;
   assign avm_channel_id_chan_Conf2Intere_active_read_data = avs_channel_id_chan_Conf2Intere_active_read_data;
   // INST channel_id_chan_Conf2Intere_active_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(1),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intere_active_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intere_active_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intere_active_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intere_active_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intere_active_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intere_active_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intere_active_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intere_active_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intere_active_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intere_active_write_almostfull = avs_channel_id_chan_Conf2Intere_active_write_almostfull;
   assign avs_channel_id_chan_Conf2Intere_cnt_write_valid = avm_channel_id_chan_Conf2Intere_cnt_write_valid;
   assign avm_channel_id_chan_Conf2Intere_cnt_write_ready = avs_channel_id_chan_Conf2Intere_cnt_write_ready;
   assign avs_channel_id_chan_Conf2Intere_cnt_write_data = avm_channel_id_chan_Conf2Intere_cnt_write_data;
   assign avm_channel_id_chan_Conf2Intere_cnt_read_valid = avs_channel_id_chan_Conf2Intere_cnt_read_valid;
   assign avs_channel_id_chan_Conf2Intere_cnt_read_ready = avm_channel_id_chan_Conf2Intere_cnt_read_ready;
   assign avm_channel_id_chan_Conf2Intere_cnt_read_data = avs_channel_id_chan_Conf2Intere_cnt_read_data;
   // INST channel_id_chan_Conf2Intere_cnt_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(1),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intere_cnt_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intere_cnt_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intere_cnt_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intere_cnt_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intere_cnt_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intere_cnt_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intere_cnt_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intere_cnt_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intere_cnt_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intere_cnt_write_almostfull = avs_channel_id_chan_Conf2Intere_cnt_write_almostfull;
   assign avs_channel_id_chan_Conf2Intere_mode_write_valid = avm_channel_id_chan_Conf2Intere_mode_write_valid;
   assign avm_channel_id_chan_Conf2Intere_mode_write_ready = avs_channel_id_chan_Conf2Intere_mode_write_ready;
   assign avs_channel_id_chan_Conf2Intere_mode_write_data = avm_channel_id_chan_Conf2Intere_mode_write_data;
   assign avm_channel_id_chan_Conf2Intere_mode_read_valid = avs_channel_id_chan_Conf2Intere_mode_read_valid;
   assign avs_channel_id_chan_Conf2Intere_mode_read_ready = avm_channel_id_chan_Conf2Intere_mode_read_ready;
   assign avm_channel_id_chan_Conf2Intere_mode_read_data = avs_channel_id_chan_Conf2Intere_mode_read_data;
   // INST channel_id_chan_Conf2Intere_mode_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(1),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intere_mode_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intere_mode_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intere_mode_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intere_mode_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intere_mode_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intere_mode_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intere_mode_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intere_mode_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intere_mode_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intere_mode_write_almostfull = avs_channel_id_chan_Conf2Intere_mode_write_almostfull;
   assign avs_channel_id_chan_Conf2Intere_x_write_valid = avm_channel_id_chan_Conf2Intere_x_write_valid;
   assign avm_channel_id_chan_Conf2Intere_x_write_ready = avs_channel_id_chan_Conf2Intere_x_write_ready;
   assign avs_channel_id_chan_Conf2Intere_x_write_data = avm_channel_id_chan_Conf2Intere_x_write_data;
   assign avm_channel_id_chan_Conf2Intere_x_read_valid = avs_channel_id_chan_Conf2Intere_x_read_valid;
   assign avs_channel_id_chan_Conf2Intere_x_read_ready = avm_channel_id_chan_Conf2Intere_x_read_ready;
   assign avm_channel_id_chan_Conf2Intere_x_read_data = avs_channel_id_chan_Conf2Intere_x_read_data;
   // INST channel_id_chan_Conf2Intere_x_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(110),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intere_x_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intere_x_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intere_x_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intere_x_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intere_x_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intere_x_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intere_x_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intere_x_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intere_x_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intere_x_write_almostfull = avs_channel_id_chan_Conf2Intere_x_write_almostfull;
   assign avs_channel_id_chan_Conf2Intere_y_write_valid = avm_channel_id_chan_Conf2Intere_y_write_valid;
   assign avm_channel_id_chan_Conf2Intere_y_write_ready = avs_channel_id_chan_Conf2Intere_y_write_ready;
   assign avs_channel_id_chan_Conf2Intere_y_write_data = avm_channel_id_chan_Conf2Intere_y_write_data;
   assign avm_channel_id_chan_Conf2Intere_y_read_valid = avs_channel_id_chan_Conf2Intere_y_read_valid;
   assign avs_channel_id_chan_Conf2Intere_y_read_ready = avm_channel_id_chan_Conf2Intere_y_read_ready;
   assign avm_channel_id_chan_Conf2Intere_y_read_data = avs_channel_id_chan_Conf2Intere_y_read_data;
   // INST channel_id_chan_Conf2Intere_y_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(100),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intere_y_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intere_y_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intere_y_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intere_y_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intere_y_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intere_y_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intere_y_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intere_y_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intere_y_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intere_y_write_almostfull = avs_channel_id_chan_Conf2Intere_y_write_almostfull;
   assign avs_channel_id_chan_Conf2Intere_z_write_valid = avm_channel_id_chan_Conf2Intere_z_write_valid;
   assign avm_channel_id_chan_Conf2Intere_z_write_ready = avs_channel_id_chan_Conf2Intere_z_write_ready;
   assign avs_channel_id_chan_Conf2Intere_z_write_data = avm_channel_id_chan_Conf2Intere_z_write_data;
   assign avm_channel_id_chan_Conf2Intere_z_read_valid = avs_channel_id_chan_Conf2Intere_z_read_valid;
   assign avs_channel_id_chan_Conf2Intere_z_read_ready = avm_channel_id_chan_Conf2Intere_z_read_ready;
   assign avm_channel_id_chan_Conf2Intere_z_read_data = avs_channel_id_chan_Conf2Intere_z_read_data;
   // INST channel_id_chan_Conf2Intere_z_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(90),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intere_z_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intere_z_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intere_z_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intere_z_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intere_z_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intere_z_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intere_z_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intere_z_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intere_z_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intere_z_write_almostfull = avs_channel_id_chan_Conf2Intere_z_write_almostfull;
   assign avs_channel_id_chan_Conf2Intrae_active_write_valid = avm_channel_id_chan_Conf2Intrae_active_write_valid;
   assign avm_channel_id_chan_Conf2Intrae_active_write_ready = avs_channel_id_chan_Conf2Intrae_active_write_ready;
   assign avs_channel_id_chan_Conf2Intrae_active_write_data = avm_channel_id_chan_Conf2Intrae_active_write_data;
   assign avm_channel_id_chan_Conf2Intrae_active_read_valid = avs_channel_id_chan_Conf2Intrae_active_read_valid;
   assign avs_channel_id_chan_Conf2Intrae_active_read_ready = avm_channel_id_chan_Conf2Intrae_active_read_ready;
   assign avm_channel_id_chan_Conf2Intrae_active_read_data = avs_channel_id_chan_Conf2Intrae_active_read_data;
   // INST channel_id_chan_Conf2Intrae_active_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(1),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intrae_active_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intrae_active_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intrae_active_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intrae_active_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intrae_active_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intrae_active_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intrae_active_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intrae_active_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intrae_active_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intrae_active_write_almostfull = avs_channel_id_chan_Conf2Intrae_active_write_almostfull;
   assign avs_channel_id_chan_Conf2Intrae_cnt_write_valid = avm_channel_id_chan_Conf2Intrae_cnt_write_valid;
   assign avm_channel_id_chan_Conf2Intrae_cnt_write_ready = avs_channel_id_chan_Conf2Intrae_cnt_write_ready;
   assign avs_channel_id_chan_Conf2Intrae_cnt_write_data = avm_channel_id_chan_Conf2Intrae_cnt_write_data;
   assign avm_channel_id_chan_Conf2Intrae_cnt_read_valid = avs_channel_id_chan_Conf2Intrae_cnt_read_valid;
   assign avs_channel_id_chan_Conf2Intrae_cnt_read_ready = avm_channel_id_chan_Conf2Intrae_cnt_read_ready;
   assign avm_channel_id_chan_Conf2Intrae_cnt_read_data = avs_channel_id_chan_Conf2Intrae_cnt_read_data;
   // INST channel_id_chan_Conf2Intrae_cnt_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(1),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intrae_cnt_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intrae_cnt_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intrae_cnt_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intrae_cnt_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intrae_cnt_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intrae_cnt_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intrae_cnt_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intrae_cnt_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intrae_cnt_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intrae_cnt_write_almostfull = avs_channel_id_chan_Conf2Intrae_cnt_write_almostfull;
   assign avs_channel_id_chan_Conf2Intrae_mode_write_valid = avm_channel_id_chan_Conf2Intrae_mode_write_valid;
   assign avm_channel_id_chan_Conf2Intrae_mode_write_ready = avs_channel_id_chan_Conf2Intrae_mode_write_ready;
   assign avs_channel_id_chan_Conf2Intrae_mode_write_data = avm_channel_id_chan_Conf2Intrae_mode_write_data;
   assign avm_channel_id_chan_Conf2Intrae_mode_read_valid = avs_channel_id_chan_Conf2Intrae_mode_read_valid;
   assign avs_channel_id_chan_Conf2Intrae_mode_read_ready = avm_channel_id_chan_Conf2Intrae_mode_read_ready;
   assign avm_channel_id_chan_Conf2Intrae_mode_read_data = avs_channel_id_chan_Conf2Intrae_mode_read_data;
   // INST channel_id_chan_Conf2Intrae_mode_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(1),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intrae_mode_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intrae_mode_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intrae_mode_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intrae_mode_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intrae_mode_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intrae_mode_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intrae_mode_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intrae_mode_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intrae_mode_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intrae_mode_write_almostfull = avs_channel_id_chan_Conf2Intrae_mode_write_almostfull;
   assign avs_channel_id_chan_Conf2Intrae_x_write_valid = avm_channel_id_chan_Conf2Intrae_x_write_valid;
   assign avm_channel_id_chan_Conf2Intrae_x_write_ready = avs_channel_id_chan_Conf2Intrae_x_write_ready;
   assign avs_channel_id_chan_Conf2Intrae_x_write_data = avm_channel_id_chan_Conf2Intrae_x_write_data;
   assign avm_channel_id_chan_Conf2Intrae_x_read_valid = avs_channel_id_chan_Conf2Intrae_x_read_valid;
   assign avs_channel_id_chan_Conf2Intrae_x_read_ready = avm_channel_id_chan_Conf2Intrae_x_read_ready;
   assign avm_channel_id_chan_Conf2Intrae_x_read_data = avs_channel_id_chan_Conf2Intrae_x_read_data;
   // INST channel_id_chan_Conf2Intrae_x_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(110),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intrae_x_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intrae_x_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intrae_x_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intrae_x_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intrae_x_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intrae_x_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intrae_x_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intrae_x_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intrae_x_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intrae_x_write_almostfull = avs_channel_id_chan_Conf2Intrae_x_write_almostfull;
   assign avs_channel_id_chan_Conf2Intrae_y_write_valid = avm_channel_id_chan_Conf2Intrae_y_write_valid;
   assign avm_channel_id_chan_Conf2Intrae_y_write_ready = avs_channel_id_chan_Conf2Intrae_y_write_ready;
   assign avs_channel_id_chan_Conf2Intrae_y_write_data = avm_channel_id_chan_Conf2Intrae_y_write_data;
   assign avm_channel_id_chan_Conf2Intrae_y_read_valid = avs_channel_id_chan_Conf2Intrae_y_read_valid;
   assign avs_channel_id_chan_Conf2Intrae_y_read_ready = avm_channel_id_chan_Conf2Intrae_y_read_ready;
   assign avm_channel_id_chan_Conf2Intrae_y_read_data = avs_channel_id_chan_Conf2Intrae_y_read_data;
   // INST channel_id_chan_Conf2Intrae_y_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(100),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intrae_y_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intrae_y_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intrae_y_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intrae_y_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intrae_y_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intrae_y_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intrae_y_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intrae_y_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intrae_y_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intrae_y_write_almostfull = avs_channel_id_chan_Conf2Intrae_y_write_almostfull;
   assign avs_channel_id_chan_Conf2Intrae_z_write_valid = avm_channel_id_chan_Conf2Intrae_z_write_valid;
   assign avm_channel_id_chan_Conf2Intrae_z_write_ready = avs_channel_id_chan_Conf2Intrae_z_write_ready;
   assign avs_channel_id_chan_Conf2Intrae_z_write_data = avm_channel_id_chan_Conf2Intrae_z_write_data;
   assign avm_channel_id_chan_Conf2Intrae_z_read_valid = avs_channel_id_chan_Conf2Intrae_z_read_valid;
   assign avs_channel_id_chan_Conf2Intrae_z_read_ready = avm_channel_id_chan_Conf2Intrae_z_read_ready;
   assign avm_channel_id_chan_Conf2Intrae_z_read_data = avs_channel_id_chan_Conf2Intrae_z_read_data;
   // INST channel_id_chan_Conf2Intrae_z_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(90),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Conf2Intrae_z_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Conf2Intrae_z_write_valid),
      .avst_in_ready(avs_channel_id_chan_Conf2Intrae_z_write_ready),
      .avst_in_data(avs_channel_id_chan_Conf2Intrae_z_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Conf2Intrae_z_read_valid),
      .avst_out_ready(avs_channel_id_chan_Conf2Intrae_z_read_ready),
      .avst_out_data(avs_channel_id_chan_Conf2Intrae_z_read_data),
      .profile_fifosize(avs_channel_id_chan_Conf2Intrae_z_fifosize),
      .almost_full(avs_channel_id_chan_Conf2Intrae_z_write_almostfull)
   );

   assign avm_channel_id_chan_Conf2Intrae_z_write_almostfull = avs_channel_id_chan_Conf2Intrae_z_write_almostfull;
   assign avs_channel_id_chan_Intere2Store_active_write_valid = avm_channel_id_chan_Intere2Store_active_write_valid;
   assign avm_channel_id_chan_Intere2Store_active_write_ready = avs_channel_id_chan_Intere2Store_active_write_ready;
   assign avs_channel_id_chan_Intere2Store_active_write_data = avm_channel_id_chan_Intere2Store_active_write_data;
   assign avm_channel_id_chan_Intere2Store_active_read_valid = avs_channel_id_chan_Intere2Store_active_read_valid;
   assign avs_channel_id_chan_Intere2Store_active_read_ready = avm_channel_id_chan_Intere2Store_active_read_ready;
   assign avm_channel_id_chan_Intere2Store_active_read_data = avs_channel_id_chan_Intere2Store_active_read_data;
   // INST channel_id_chan_Intere2Store_active_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(392),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intere2Store_active_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intere2Store_active_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intere2Store_active_write_ready),
      .avst_in_data(avs_channel_id_chan_Intere2Store_active_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intere2Store_active_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intere2Store_active_read_ready),
      .avst_out_data(avs_channel_id_chan_Intere2Store_active_read_data),
      .profile_fifosize(avs_channel_id_chan_Intere2Store_active_fifosize),
      .almost_full(avs_channel_id_chan_Intere2Store_active_write_almostfull)
   );

   assign avm_channel_id_chan_Intere2Store_active_write_almostfull = avs_channel_id_chan_Intere2Store_active_write_almostfull;
   assign avs_channel_id_chan_Intere2Store_cnt_write_valid = avm_channel_id_chan_Intere2Store_cnt_write_valid;
   assign avm_channel_id_chan_Intere2Store_cnt_write_ready = avs_channel_id_chan_Intere2Store_cnt_write_ready;
   assign avs_channel_id_chan_Intere2Store_cnt_write_data = avm_channel_id_chan_Intere2Store_cnt_write_data;
   assign avm_channel_id_chan_Intere2Store_cnt_read_valid = avs_channel_id_chan_Intere2Store_cnt_read_valid;
   assign avs_channel_id_chan_Intere2Store_cnt_read_ready = avm_channel_id_chan_Intere2Store_cnt_read_ready;
   assign avm_channel_id_chan_Intere2Store_cnt_read_data = avs_channel_id_chan_Intere2Store_cnt_read_data;
   // INST channel_id_chan_Intere2Store_cnt_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(391),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intere2Store_cnt_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intere2Store_cnt_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intere2Store_cnt_write_ready),
      .avst_in_data(avs_channel_id_chan_Intere2Store_cnt_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intere2Store_cnt_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intere2Store_cnt_read_ready),
      .avst_out_data(avs_channel_id_chan_Intere2Store_cnt_read_data),
      .profile_fifosize(avs_channel_id_chan_Intere2Store_cnt_fifosize),
      .almost_full(avs_channel_id_chan_Intere2Store_cnt_write_almostfull)
   );

   assign avm_channel_id_chan_Intere2Store_cnt_write_almostfull = avs_channel_id_chan_Intere2Store_cnt_write_almostfull;
   assign avs_channel_id_chan_Intere2Store_intere_write_valid = avm_channel_id_chan_Intere2Store_intere_write_valid;
   assign avm_channel_id_chan_Intere2Store_intere_write_ready = avs_channel_id_chan_Intere2Store_intere_write_ready;
   assign avs_channel_id_chan_Intere2Store_intere_write_data = avm_channel_id_chan_Intere2Store_intere_write_data;
   assign avm_channel_id_chan_Intere2Store_intere_read_valid = avs_channel_id_chan_Intere2Store_intere_read_valid;
   assign avs_channel_id_chan_Intere2Store_intere_read_ready = avm_channel_id_chan_Intere2Store_intere_read_ready;
   assign avm_channel_id_chan_Intere2Store_intere_read_data = avs_channel_id_chan_Intere2Store_intere_read_data;
   // INST channel_id_chan_Intere2Store_intere_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(393),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intere2Store_intere_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intere2Store_intere_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intere2Store_intere_write_ready),
      .avst_in_data(avs_channel_id_chan_Intere2Store_intere_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intere2Store_intere_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intere2Store_intere_read_ready),
      .avst_out_data(avs_channel_id_chan_Intere2Store_intere_read_data),
      .profile_fifosize(avs_channel_id_chan_Intere2Store_intere_fifosize),
      .almost_full(avs_channel_id_chan_Intere2Store_intere_write_almostfull)
   );

   assign avm_channel_id_chan_Intere2Store_intere_write_almostfull = avs_channel_id_chan_Intere2Store_intere_write_almostfull;
   assign avs_channel_id_chan_Intere2Store_mode_write_valid = avm_channel_id_chan_Intere2Store_mode_write_valid;
   assign avm_channel_id_chan_Intere2Store_mode_write_ready = avs_channel_id_chan_Intere2Store_mode_write_ready;
   assign avs_channel_id_chan_Intere2Store_mode_write_data = avm_channel_id_chan_Intere2Store_mode_write_data;
   assign avm_channel_id_chan_Intere2Store_mode_read_valid = avs_channel_id_chan_Intere2Store_mode_read_valid;
   assign avs_channel_id_chan_Intere2Store_mode_read_ready = avm_channel_id_chan_Intere2Store_mode_read_ready;
   assign avm_channel_id_chan_Intere2Store_mode_read_data = avs_channel_id_chan_Intere2Store_mode_read_data;
   // INST channel_id_chan_Intere2Store_mode_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(391),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intere2Store_mode_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intere2Store_mode_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intere2Store_mode_write_ready),
      .avst_in_data(avs_channel_id_chan_Intere2Store_mode_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intere2Store_mode_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intere2Store_mode_read_ready),
      .avst_out_data(avs_channel_id_chan_Intere2Store_mode_read_data),
      .profile_fifosize(avs_channel_id_chan_Intere2Store_mode_fifosize),
      .almost_full(avs_channel_id_chan_Intere2Store_mode_write_almostfull)
   );

   assign avm_channel_id_chan_Intere2Store_mode_write_almostfull = avs_channel_id_chan_Intere2Store_mode_write_almostfull;
   assign avs_channel_id_chan_Intrae2Store_active_write_valid = avm_channel_id_chan_Intrae2Store_active_write_valid;
   assign avm_channel_id_chan_Intrae2Store_active_write_ready = avs_channel_id_chan_Intrae2Store_active_write_ready;
   assign avs_channel_id_chan_Intrae2Store_active_write_data = avm_channel_id_chan_Intrae2Store_active_write_data;
   assign avm_channel_id_chan_Intrae2Store_active_read_valid = avs_channel_id_chan_Intrae2Store_active_read_valid;
   assign avs_channel_id_chan_Intrae2Store_active_read_ready = avm_channel_id_chan_Intrae2Store_active_read_ready;
   assign avm_channel_id_chan_Intrae2Store_active_read_data = avs_channel_id_chan_Intrae2Store_active_read_data;
   // INST channel_id_chan_Intrae2Store_active_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(0),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intrae2Store_active_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intrae2Store_active_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intrae2Store_active_write_ready),
      .avst_in_data(avs_channel_id_chan_Intrae2Store_active_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intrae2Store_active_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intrae2Store_active_read_ready),
      .avst_out_data(avs_channel_id_chan_Intrae2Store_active_read_data),
      .profile_fifosize(avs_channel_id_chan_Intrae2Store_active_fifosize),
      .almost_full(avs_channel_id_chan_Intrae2Store_active_write_almostfull)
   );

   assign avm_channel_id_chan_Intrae2Store_active_write_almostfull = avs_channel_id_chan_Intrae2Store_active_write_almostfull;
   assign avs_channel_id_chan_Intrae2Store_cnt_write_valid = avm_channel_id_chan_Intrae2Store_cnt_write_valid;
   assign avm_channel_id_chan_Intrae2Store_cnt_write_ready = avs_channel_id_chan_Intrae2Store_cnt_write_ready;
   assign avs_channel_id_chan_Intrae2Store_cnt_write_data = avm_channel_id_chan_Intrae2Store_cnt_write_data;
   assign avm_channel_id_chan_Intrae2Store_cnt_read_valid = avs_channel_id_chan_Intrae2Store_cnt_read_valid;
   assign avs_channel_id_chan_Intrae2Store_cnt_read_ready = avm_channel_id_chan_Intrae2Store_cnt_read_ready;
   assign avm_channel_id_chan_Intrae2Store_cnt_read_data = avs_channel_id_chan_Intrae2Store_cnt_read_data;
   // INST channel_id_chan_Intrae2Store_cnt_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(0),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intrae2Store_cnt_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intrae2Store_cnt_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intrae2Store_cnt_write_ready),
      .avst_in_data(avs_channel_id_chan_Intrae2Store_cnt_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intrae2Store_cnt_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intrae2Store_cnt_read_ready),
      .avst_out_data(avs_channel_id_chan_Intrae2Store_cnt_read_data),
      .profile_fifosize(avs_channel_id_chan_Intrae2Store_cnt_fifosize),
      .almost_full(avs_channel_id_chan_Intrae2Store_cnt_write_almostfull)
   );

   assign avm_channel_id_chan_Intrae2Store_cnt_write_almostfull = avs_channel_id_chan_Intrae2Store_cnt_write_almostfull;
   assign avs_channel_id_chan_Intrae2Store_intrae_write_valid = avm_channel_id_chan_Intrae2Store_intrae_write_valid;
   assign avm_channel_id_chan_Intrae2Store_intrae_write_ready = avs_channel_id_chan_Intrae2Store_intrae_write_ready;
   assign avs_channel_id_chan_Intrae2Store_intrae_write_data = avm_channel_id_chan_Intrae2Store_intrae_write_data;
   assign avm_channel_id_chan_Intrae2Store_intrae_read_valid = avs_channel_id_chan_Intrae2Store_intrae_read_valid;
   assign avs_channel_id_chan_Intrae2Store_intrae_read_ready = avm_channel_id_chan_Intrae2Store_intrae_read_ready;
   assign avm_channel_id_chan_Intrae2Store_intrae_read_data = avs_channel_id_chan_Intrae2Store_intrae_read_data;
   // INST channel_id_chan_Intrae2Store_intrae_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(1),
      .DATA_W(32),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intrae2Store_intrae_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intrae2Store_intrae_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intrae2Store_intrae_write_ready),
      .avst_in_data(avs_channel_id_chan_Intrae2Store_intrae_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intrae2Store_intrae_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intrae2Store_intrae_read_ready),
      .avst_out_data(avs_channel_id_chan_Intrae2Store_intrae_read_data),
      .profile_fifosize(avs_channel_id_chan_Intrae2Store_intrae_fifosize),
      .almost_full(avs_channel_id_chan_Intrae2Store_intrae_write_almostfull)
   );

   assign avm_channel_id_chan_Intrae2Store_intrae_write_almostfull = avs_channel_id_chan_Intrae2Store_intrae_write_almostfull;
   assign avs_channel_id_chan_Intrae2Store_mode_write_valid = avm_channel_id_chan_Intrae2Store_mode_write_valid;
   assign avm_channel_id_chan_Intrae2Store_mode_write_ready = avs_channel_id_chan_Intrae2Store_mode_write_ready;
   assign avs_channel_id_chan_Intrae2Store_mode_write_data = avm_channel_id_chan_Intrae2Store_mode_write_data;
   assign avm_channel_id_chan_Intrae2Store_mode_read_valid = avs_channel_id_chan_Intrae2Store_mode_read_valid;
   assign avs_channel_id_chan_Intrae2Store_mode_read_ready = avm_channel_id_chan_Intrae2Store_mode_read_ready;
   assign avm_channel_id_chan_Intrae2Store_mode_read_data = avs_channel_id_chan_Intrae2Store_mode_read_data;
   // INST channel_id_chan_Intrae2Store_mode_fifo of acl_channel_fifo
   acl_channel_fifo
   #(
      .FIFO_DEPTH(0),
      .DATA_W(8),
      .ADJUST_FOR_LATENCY(1),
      .INTENDED_DEVICE_FAMILY("Arria 10")
   )
   channel_id_chan_Intrae2Store_mode_fifo
   (
      .clock(clock),
      .resetn(resetn),
      // AVST avst_in
      .avst_in_valid(avs_channel_id_chan_Intrae2Store_mode_write_valid),
      .avst_in_ready(avs_channel_id_chan_Intrae2Store_mode_write_ready),
      .avst_in_data(avs_channel_id_chan_Intrae2Store_mode_write_data),
      // AVST avst_out
      .avst_out_valid(avs_channel_id_chan_Intrae2Store_mode_read_valid),
      .avst_out_ready(avs_channel_id_chan_Intrae2Store_mode_read_ready),
      .avst_out_data(avs_channel_id_chan_Intrae2Store_mode_read_data),
      .profile_fifosize(avs_channel_id_chan_Intrae2Store_mode_fifosize),
      .almost_full(avs_channel_id_chan_Intrae2Store_mode_write_almostfull)
   );

   assign avm_channel_id_chan_Intrae2Store_mode_write_almostfull = avs_channel_id_chan_Intrae2Store_mode_write_almostfull;
   assign counter_resetn_Krnl_Conform = (resetn & ~(Krnl_Conform_start));
   // INST Krnl_Conform_printf_counter of acl_printf_counter
   acl_printf_counter Krnl_Conform_printf_counter
   (
      .clock(clock),
      .resetn(counter_resetn_Krnl_Conform),
      .enable(counter_enable_Krnl_Conform),
      .i_counter_reset(counter_reset_Krnl_Conform),
      .i_increment(counter_increment_Krnl_Conform),
      .i_init(counter_init_Krnl_Conform),
      .i_limit(counter_limit_Krnl_Conform),
      .o_size(counter_size_Krnl_Conform),
      .o_resultvalid(counter_resultvalid_Krnl_Conform),
      .o_result(counter_result_Krnl_Conform),
      .o_full(counter_full_Krnl_Conform),
      .o_stall(counter_stall_Krnl_Conform)
   );

   generate
   begin:printf_system_Krnl_Conform
      logic p_icm_arb_in_arb_request [1];
      logic p_icm_arb_in_arb_enable [1];
      logic p_icm_arb_in_arb_read [1];
      logic p_icm_arb_in_arb_write [1];
      logic [5:0] p_icm_arb_in_arb_burstcount [1];
      logic [26:0] p_icm_arb_in_arb_address [1];
      logic [255:0] p_icm_arb_in_arb_writedata [1];
      logic [31:0] p_icm_arb_in_arb_byteenable [1];
      logic p_icm_arb_in_arb_stall [1];
      logic p_icm_arb_in_wrp_ack [1];
      logic p_icm_arb_in_rrp_datavalid [1];
      logic [255:0] p_icm_arb_in_rrp_data [1];
      logic p_icm_arb_out_arb_request [1];
      logic p_icm_arb_out_arb_enable [1];
      logic p_icm_arb_out_arb_read [1];
      logic p_icm_arb_out_arb_write [1];
      logic [5:0] p_icm_arb_out_arb_burstcount [1];
      logic [26:0] p_icm_arb_out_arb_address [1];
      logic [255:0] p_icm_arb_out_arb_writedata [1];
      logic [31:0] p_icm_arb_out_arb_byteenable [1];
      logic p_icm_arb_out_arb_stall [1];
      logic p_icm_arb_out_wrp_ack [1];
      logic p_icm_arb_out_rrp_datavalid [1];
      logic [255:0] p_icm_arb_out_rrp_data [1];

      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:p
         // INST p_avm_to_ic of acl_avm_to_ic
         acl_avm_to_ic
         #(
            .DATA_W(256),
            .WRITEDATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(32),
            .BYTEENA_W(32),
            .ADDR_SHIFT(0)
         )
         p_avm_to_ic
         (
            // AVM avm
            .avm_enable(printf_avm_Krnl_Conform_enable[__i]),
            .avm_read(printf_avm_Krnl_Conform_read[__i]),
            .avm_write(printf_avm_Krnl_Conform_write[__i]),
            .avm_burstcount(printf_avm_Krnl_Conform_burstcount[__i]),
            .avm_address(printf_avm_Krnl_Conform_address[__i]),
            .avm_writedata(printf_avm_Krnl_Conform_writedata[__i]),
            .avm_byteenable(printf_avm_Krnl_Conform_byteenable[__i]),
            .avm_waitrequest(printf_avm_Krnl_Conform_waitrequest[__i]),
            .avm_readdata(printf_avm_Krnl_Conform_readdata[__i]),
            .avm_readdatavalid(printf_avm_Krnl_Conform_readdatavalid[__i]),
            // ICM ic
            .ic_arb_request(p_icm_arb_in_arb_request[__i]),
            .ic_arb_enable(p_icm_arb_in_arb_enable[__i]),
            .ic_arb_read(p_icm_arb_in_arb_read[__i]),
            .ic_arb_write(p_icm_arb_in_arb_write[__i]),
            .ic_arb_burstcount(p_icm_arb_in_arb_burstcount[__i]),
            .ic_arb_address(p_icm_arb_in_arb_address[__i]),
            .ic_arb_writedata(p_icm_arb_in_arb_writedata[__i]),
            .ic_arb_byteenable(p_icm_arb_in_arb_byteenable[__i]),
            .ic_arb_stall(p_icm_arb_in_arb_stall[__i]),
            .ic_wrp_ack(p_icm_arb_in_wrp_ack[__i]),
            .ic_rrp_datavalid(p_icm_arb_in_rrp_datavalid[__i]),
            .ic_rrp_data(p_icm_arb_in_rrp_data[__i])
         );

      end

      // INST printf_ic_interconnect of Krnl_GA_system_interconnect_6
      Krnl_GA_system_interconnect_6 printf_ic_interconnect
      (
         .clock(clock),
         .resetn(resetn),
         // ICM m
         .m_arb_request(p_icm_arb_in_arb_request),
         .m_arb_enable(p_icm_arb_in_arb_enable),
         .m_arb_read(p_icm_arb_in_arb_read),
         .m_arb_write(p_icm_arb_in_arb_write),
         .m_arb_burstcount(p_icm_arb_in_arb_burstcount),
         .m_arb_address(p_icm_arb_in_arb_address),
         .m_arb_writedata(p_icm_arb_in_arb_writedata),
         .m_arb_byteenable(p_icm_arb_in_arb_byteenable),
         .m_arb_stall(p_icm_arb_in_arb_stall),
         .m_wrp_ack(p_icm_arb_in_wrp_ack),
         .m_rrp_datavalid(p_icm_arb_in_rrp_datavalid),
         .m_rrp_data(p_icm_arb_in_rrp_data),
         // ICM mout
         .mout_arb_request(p_icm_arb_out_arb_request[0]),
         .mout_arb_enable(p_icm_arb_out_arb_enable[0]),
         .mout_arb_read(p_icm_arb_out_arb_read[0]),
         .mout_arb_write(p_icm_arb_out_arb_write[0]),
         .mout_arb_burstcount(p_icm_arb_out_arb_burstcount[0]),
         .mout_arb_address(p_icm_arb_out_arb_address[0]),
         .mout_arb_writedata(p_icm_arb_out_arb_writedata[0]),
         .mout_arb_byteenable(p_icm_arb_out_arb_byteenable[0]),
         .mout_arb_id(),
         .mout_arb_stall(p_icm_arb_out_arb_stall[0]),
         .mout_wrp_ack(p_icm_arb_out_wrp_ack[0]),
         .mout_rrp_datavalid(p_icm_arb_out_rrp_datavalid[0]),
         .mout_rrp_data(p_icm_arb_out_rrp_data[0])
      );

      assign counter_increment_Krnl_Conform = p_icm_arb_out_arb_address[0];
      assign p_icm_arb_out_rrp_data[0] = counter_result_Krnl_Conform;
      assign counter_enable_Krnl_Conform = p_icm_arb_out_arb_read[0];
      assign p_icm_arb_out_rrp_datavalid[0] = counter_resultvalid_Krnl_Conform;
      assign p_icm_arb_out_arb_stall[0] = counter_stall_Krnl_Conform;
      assign p_icm_arb_out_wrp_ack[0] = 1'b1;
   end
   endgenerate

   assign counter_resetn_Krnl_GA = (resetn & ~(Krnl_GA_start));
   // INST Krnl_GA_printf_counter of acl_printf_counter
   acl_printf_counter Krnl_GA_printf_counter
   (
      .clock(clock),
      .resetn(counter_resetn_Krnl_GA),
      .enable(counter_enable_Krnl_GA),
      .i_counter_reset(counter_reset_Krnl_GA),
      .i_increment(counter_increment_Krnl_GA),
      .i_init(counter_init_Krnl_GA),
      .i_limit(counter_limit_Krnl_GA),
      .o_size(counter_size_Krnl_GA),
      .o_resultvalid(counter_resultvalid_Krnl_GA),
      .o_result(counter_result_Krnl_GA),
      .o_full(counter_full_Krnl_GA),
      .o_stall(counter_stall_Krnl_GA)
   );

   generate
   begin:printf_system_Krnl_GA
      logic p_icm_arb_in_arb_request [1];
      logic p_icm_arb_in_arb_enable [1];
      logic p_icm_arb_in_arb_read [1];
      logic p_icm_arb_in_arb_write [1];
      logic [5:0] p_icm_arb_in_arb_burstcount [1];
      logic [26:0] p_icm_arb_in_arb_address [1];
      logic [255:0] p_icm_arb_in_arb_writedata [1];
      logic [31:0] p_icm_arb_in_arb_byteenable [1];
      logic p_icm_arb_in_arb_stall [1];
      logic p_icm_arb_in_wrp_ack [1];
      logic p_icm_arb_in_rrp_datavalid [1];
      logic [255:0] p_icm_arb_in_rrp_data [1];
      logic p_icm_arb_out_arb_request [1];
      logic p_icm_arb_out_arb_enable [1];
      logic p_icm_arb_out_arb_read [1];
      logic p_icm_arb_out_arb_write [1];
      logic [5:0] p_icm_arb_out_arb_burstcount [1];
      logic [26:0] p_icm_arb_out_arb_address [1];
      logic [255:0] p_icm_arb_out_arb_writedata [1];
      logic [31:0] p_icm_arb_out_arb_byteenable [1];
      logic p_icm_arb_out_arb_stall [1];
      logic p_icm_arb_out_wrp_ack [1];
      logic p_icm_arb_out_rrp_datavalid [1];
      logic [255:0] p_icm_arb_out_rrp_data [1];

      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:p
         // INST p_avm_to_ic of acl_avm_to_ic
         acl_avm_to_ic
         #(
            .DATA_W(256),
            .WRITEDATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(32),
            .BYTEENA_W(32),
            .ADDR_SHIFT(0)
         )
         p_avm_to_ic
         (
            // AVM avm
            .avm_enable(printf_avm_Krnl_GA_enable[__i]),
            .avm_read(printf_avm_Krnl_GA_read[__i]),
            .avm_write(printf_avm_Krnl_GA_write[__i]),
            .avm_burstcount(printf_avm_Krnl_GA_burstcount[__i]),
            .avm_address(printf_avm_Krnl_GA_address[__i]),
            .avm_writedata(printf_avm_Krnl_GA_writedata[__i]),
            .avm_byteenable(printf_avm_Krnl_GA_byteenable[__i]),
            .avm_waitrequest(printf_avm_Krnl_GA_waitrequest[__i]),
            .avm_readdata(printf_avm_Krnl_GA_readdata[__i]),
            .avm_readdatavalid(printf_avm_Krnl_GA_readdatavalid[__i]),
            // ICM ic
            .ic_arb_request(p_icm_arb_in_arb_request[__i]),
            .ic_arb_enable(p_icm_arb_in_arb_enable[__i]),
            .ic_arb_read(p_icm_arb_in_arb_read[__i]),
            .ic_arb_write(p_icm_arb_in_arb_write[__i]),
            .ic_arb_burstcount(p_icm_arb_in_arb_burstcount[__i]),
            .ic_arb_address(p_icm_arb_in_arb_address[__i]),
            .ic_arb_writedata(p_icm_arb_in_arb_writedata[__i]),
            .ic_arb_byteenable(p_icm_arb_in_arb_byteenable[__i]),
            .ic_arb_stall(p_icm_arb_in_arb_stall[__i]),
            .ic_wrp_ack(p_icm_arb_in_wrp_ack[__i]),
            .ic_rrp_datavalid(p_icm_arb_in_rrp_datavalid[__i]),
            .ic_rrp_data(p_icm_arb_in_rrp_data[__i])
         );

      end

      // INST printf_ic_interconnect of Krnl_GA_system_interconnect_6
      Krnl_GA_system_interconnect_6 printf_ic_interconnect
      (
         .clock(clock),
         .resetn(resetn),
         // ICM m
         .m_arb_request(p_icm_arb_in_arb_request),
         .m_arb_enable(p_icm_arb_in_arb_enable),
         .m_arb_read(p_icm_arb_in_arb_read),
         .m_arb_write(p_icm_arb_in_arb_write),
         .m_arb_burstcount(p_icm_arb_in_arb_burstcount),
         .m_arb_address(p_icm_arb_in_arb_address),
         .m_arb_writedata(p_icm_arb_in_arb_writedata),
         .m_arb_byteenable(p_icm_arb_in_arb_byteenable),
         .m_arb_stall(p_icm_arb_in_arb_stall),
         .m_wrp_ack(p_icm_arb_in_wrp_ack),
         .m_rrp_datavalid(p_icm_arb_in_rrp_datavalid),
         .m_rrp_data(p_icm_arb_in_rrp_data),
         // ICM mout
         .mout_arb_request(p_icm_arb_out_arb_request[0]),
         .mout_arb_enable(p_icm_arb_out_arb_enable[0]),
         .mout_arb_read(p_icm_arb_out_arb_read[0]),
         .mout_arb_write(p_icm_arb_out_arb_write[0]),
         .mout_arb_burstcount(p_icm_arb_out_arb_burstcount[0]),
         .mout_arb_address(p_icm_arb_out_arb_address[0]),
         .mout_arb_writedata(p_icm_arb_out_arb_writedata[0]),
         .mout_arb_byteenable(p_icm_arb_out_arb_byteenable[0]),
         .mout_arb_id(),
         .mout_arb_stall(p_icm_arb_out_arb_stall[0]),
         .mout_wrp_ack(p_icm_arb_out_wrp_ack[0]),
         .mout_rrp_datavalid(p_icm_arb_out_rrp_datavalid[0]),
         .mout_rrp_data(p_icm_arb_out_rrp_data[0])
      );

      assign counter_increment_Krnl_GA = p_icm_arb_out_arb_address[0];
      assign p_icm_arb_out_rrp_data[0] = counter_result_Krnl_GA;
      assign counter_enable_Krnl_GA = p_icm_arb_out_arb_read[0];
      assign p_icm_arb_out_rrp_datavalid[0] = counter_resultvalid_Krnl_GA;
      assign p_icm_arb_out_arb_stall[0] = counter_stall_Krnl_GA;
      assign p_icm_arb_out_wrp_ack[0] = 1'b1;
   end
   endgenerate

   assign counter_resetn_Krnl_InterE = (resetn & ~(Krnl_InterE_start));
   // INST Krnl_InterE_printf_counter of acl_printf_counter
   acl_printf_counter Krnl_InterE_printf_counter
   (
      .clock(clock),
      .resetn(counter_resetn_Krnl_InterE),
      .enable(counter_enable_Krnl_InterE),
      .i_counter_reset(counter_reset_Krnl_InterE),
      .i_increment(counter_increment_Krnl_InterE),
      .i_init(counter_init_Krnl_InterE),
      .i_limit(counter_limit_Krnl_InterE),
      .o_size(counter_size_Krnl_InterE),
      .o_resultvalid(counter_resultvalid_Krnl_InterE),
      .o_result(counter_result_Krnl_InterE),
      .o_full(counter_full_Krnl_InterE),
      .o_stall(counter_stall_Krnl_InterE)
   );

   generate
   begin:printf_system_Krnl_InterE
      logic p_icm_arb_in_arb_request [1];
      logic p_icm_arb_in_arb_enable [1];
      logic p_icm_arb_in_arb_read [1];
      logic p_icm_arb_in_arb_write [1];
      logic [5:0] p_icm_arb_in_arb_burstcount [1];
      logic [26:0] p_icm_arb_in_arb_address [1];
      logic [255:0] p_icm_arb_in_arb_writedata [1];
      logic [31:0] p_icm_arb_in_arb_byteenable [1];
      logic p_icm_arb_in_arb_stall [1];
      logic p_icm_arb_in_wrp_ack [1];
      logic p_icm_arb_in_rrp_datavalid [1];
      logic [255:0] p_icm_arb_in_rrp_data [1];
      logic p_icm_arb_out_arb_request [1];
      logic p_icm_arb_out_arb_enable [1];
      logic p_icm_arb_out_arb_read [1];
      logic p_icm_arb_out_arb_write [1];
      logic [5:0] p_icm_arb_out_arb_burstcount [1];
      logic [26:0] p_icm_arb_out_arb_address [1];
      logic [255:0] p_icm_arb_out_arb_writedata [1];
      logic [31:0] p_icm_arb_out_arb_byteenable [1];
      logic p_icm_arb_out_arb_stall [1];
      logic p_icm_arb_out_wrp_ack [1];
      logic p_icm_arb_out_rrp_datavalid [1];
      logic [255:0] p_icm_arb_out_rrp_data [1];

      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:p
         // INST p_avm_to_ic of acl_avm_to_ic
         acl_avm_to_ic
         #(
            .DATA_W(256),
            .WRITEDATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(32),
            .BYTEENA_W(32),
            .ADDR_SHIFT(0)
         )
         p_avm_to_ic
         (
            // AVM avm
            .avm_enable(printf_avm_Krnl_InterE_enable[__i]),
            .avm_read(printf_avm_Krnl_InterE_read[__i]),
            .avm_write(printf_avm_Krnl_InterE_write[__i]),
            .avm_burstcount(printf_avm_Krnl_InterE_burstcount[__i]),
            .avm_address(printf_avm_Krnl_InterE_address[__i]),
            .avm_writedata(printf_avm_Krnl_InterE_writedata[__i]),
            .avm_byteenable(printf_avm_Krnl_InterE_byteenable[__i]),
            .avm_waitrequest(printf_avm_Krnl_InterE_waitrequest[__i]),
            .avm_readdata(printf_avm_Krnl_InterE_readdata[__i]),
            .avm_readdatavalid(printf_avm_Krnl_InterE_readdatavalid[__i]),
            // ICM ic
            .ic_arb_request(p_icm_arb_in_arb_request[__i]),
            .ic_arb_enable(p_icm_arb_in_arb_enable[__i]),
            .ic_arb_read(p_icm_arb_in_arb_read[__i]),
            .ic_arb_write(p_icm_arb_in_arb_write[__i]),
            .ic_arb_burstcount(p_icm_arb_in_arb_burstcount[__i]),
            .ic_arb_address(p_icm_arb_in_arb_address[__i]),
            .ic_arb_writedata(p_icm_arb_in_arb_writedata[__i]),
            .ic_arb_byteenable(p_icm_arb_in_arb_byteenable[__i]),
            .ic_arb_stall(p_icm_arb_in_arb_stall[__i]),
            .ic_wrp_ack(p_icm_arb_in_wrp_ack[__i]),
            .ic_rrp_datavalid(p_icm_arb_in_rrp_datavalid[__i]),
            .ic_rrp_data(p_icm_arb_in_rrp_data[__i])
         );

      end

      // INST printf_ic_interconnect of Krnl_GA_system_interconnect_6
      Krnl_GA_system_interconnect_6 printf_ic_interconnect
      (
         .clock(clock),
         .resetn(resetn),
         // ICM m
         .m_arb_request(p_icm_arb_in_arb_request),
         .m_arb_enable(p_icm_arb_in_arb_enable),
         .m_arb_read(p_icm_arb_in_arb_read),
         .m_arb_write(p_icm_arb_in_arb_write),
         .m_arb_burstcount(p_icm_arb_in_arb_burstcount),
         .m_arb_address(p_icm_arb_in_arb_address),
         .m_arb_writedata(p_icm_arb_in_arb_writedata),
         .m_arb_byteenable(p_icm_arb_in_arb_byteenable),
         .m_arb_stall(p_icm_arb_in_arb_stall),
         .m_wrp_ack(p_icm_arb_in_wrp_ack),
         .m_rrp_datavalid(p_icm_arb_in_rrp_datavalid),
         .m_rrp_data(p_icm_arb_in_rrp_data),
         // ICM mout
         .mout_arb_request(p_icm_arb_out_arb_request[0]),
         .mout_arb_enable(p_icm_arb_out_arb_enable[0]),
         .mout_arb_read(p_icm_arb_out_arb_read[0]),
         .mout_arb_write(p_icm_arb_out_arb_write[0]),
         .mout_arb_burstcount(p_icm_arb_out_arb_burstcount[0]),
         .mout_arb_address(p_icm_arb_out_arb_address[0]),
         .mout_arb_writedata(p_icm_arb_out_arb_writedata[0]),
         .mout_arb_byteenable(p_icm_arb_out_arb_byteenable[0]),
         .mout_arb_id(),
         .mout_arb_stall(p_icm_arb_out_arb_stall[0]),
         .mout_wrp_ack(p_icm_arb_out_wrp_ack[0]),
         .mout_rrp_datavalid(p_icm_arb_out_rrp_datavalid[0]),
         .mout_rrp_data(p_icm_arb_out_rrp_data[0])
      );

      assign counter_increment_Krnl_InterE = p_icm_arb_out_arb_address[0];
      assign p_icm_arb_out_rrp_data[0] = counter_result_Krnl_InterE;
      assign counter_enable_Krnl_InterE = p_icm_arb_out_arb_read[0];
      assign p_icm_arb_out_rrp_datavalid[0] = counter_resultvalid_Krnl_InterE;
      assign p_icm_arb_out_arb_stall[0] = counter_stall_Krnl_InterE;
      assign p_icm_arb_out_wrp_ack[0] = 1'b1;
   end
   endgenerate

   assign counter_resetn_Krnl_IntraE = (resetn & ~(Krnl_IntraE_start));
   // INST Krnl_IntraE_printf_counter of acl_printf_counter
   acl_printf_counter Krnl_IntraE_printf_counter
   (
      .clock(clock),
      .resetn(counter_resetn_Krnl_IntraE),
      .enable(counter_enable_Krnl_IntraE),
      .i_counter_reset(counter_reset_Krnl_IntraE),
      .i_increment(counter_increment_Krnl_IntraE),
      .i_init(counter_init_Krnl_IntraE),
      .i_limit(counter_limit_Krnl_IntraE),
      .o_size(counter_size_Krnl_IntraE),
      .o_resultvalid(counter_resultvalid_Krnl_IntraE),
      .o_result(counter_result_Krnl_IntraE),
      .o_full(counter_full_Krnl_IntraE),
      .o_stall(counter_stall_Krnl_IntraE)
   );

   generate
   begin:printf_system_Krnl_IntraE
      logic p_icm_arb_in_arb_request [1];
      logic p_icm_arb_in_arb_enable [1];
      logic p_icm_arb_in_arb_read [1];
      logic p_icm_arb_in_arb_write [1];
      logic [5:0] p_icm_arb_in_arb_burstcount [1];
      logic [26:0] p_icm_arb_in_arb_address [1];
      logic [255:0] p_icm_arb_in_arb_writedata [1];
      logic [31:0] p_icm_arb_in_arb_byteenable [1];
      logic p_icm_arb_in_arb_stall [1];
      logic p_icm_arb_in_wrp_ack [1];
      logic p_icm_arb_in_rrp_datavalid [1];
      logic [255:0] p_icm_arb_in_rrp_data [1];
      logic p_icm_arb_out_arb_request [1];
      logic p_icm_arb_out_arb_enable [1];
      logic p_icm_arb_out_arb_read [1];
      logic p_icm_arb_out_arb_write [1];
      logic [5:0] p_icm_arb_out_arb_burstcount [1];
      logic [26:0] p_icm_arb_out_arb_address [1];
      logic [255:0] p_icm_arb_out_arb_writedata [1];
      logic [31:0] p_icm_arb_out_arb_byteenable [1];
      logic p_icm_arb_out_arb_stall [1];
      logic p_icm_arb_out_wrp_ack [1];
      logic p_icm_arb_out_rrp_datavalid [1];
      logic [255:0] p_icm_arb_out_rrp_data [1];

      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:p
         // INST p_avm_to_ic of acl_avm_to_ic
         acl_avm_to_ic
         #(
            .DATA_W(256),
            .WRITEDATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(32),
            .BYTEENA_W(32),
            .ADDR_SHIFT(0)
         )
         p_avm_to_ic
         (
            // AVM avm
            .avm_enable(printf_avm_Krnl_IntraE_enable[__i]),
            .avm_read(printf_avm_Krnl_IntraE_read[__i]),
            .avm_write(printf_avm_Krnl_IntraE_write[__i]),
            .avm_burstcount(printf_avm_Krnl_IntraE_burstcount[__i]),
            .avm_address(printf_avm_Krnl_IntraE_address[__i]),
            .avm_writedata(printf_avm_Krnl_IntraE_writedata[__i]),
            .avm_byteenable(printf_avm_Krnl_IntraE_byteenable[__i]),
            .avm_waitrequest(printf_avm_Krnl_IntraE_waitrequest[__i]),
            .avm_readdata(printf_avm_Krnl_IntraE_readdata[__i]),
            .avm_readdatavalid(printf_avm_Krnl_IntraE_readdatavalid[__i]),
            // ICM ic
            .ic_arb_request(p_icm_arb_in_arb_request[__i]),
            .ic_arb_enable(p_icm_arb_in_arb_enable[__i]),
            .ic_arb_read(p_icm_arb_in_arb_read[__i]),
            .ic_arb_write(p_icm_arb_in_arb_write[__i]),
            .ic_arb_burstcount(p_icm_arb_in_arb_burstcount[__i]),
            .ic_arb_address(p_icm_arb_in_arb_address[__i]),
            .ic_arb_writedata(p_icm_arb_in_arb_writedata[__i]),
            .ic_arb_byteenable(p_icm_arb_in_arb_byteenable[__i]),
            .ic_arb_stall(p_icm_arb_in_arb_stall[__i]),
            .ic_wrp_ack(p_icm_arb_in_wrp_ack[__i]),
            .ic_rrp_datavalid(p_icm_arb_in_rrp_datavalid[__i]),
            .ic_rrp_data(p_icm_arb_in_rrp_data[__i])
         );

      end

      // INST printf_ic_interconnect of Krnl_GA_system_interconnect_6
      Krnl_GA_system_interconnect_6 printf_ic_interconnect
      (
         .clock(clock),
         .resetn(resetn),
         // ICM m
         .m_arb_request(p_icm_arb_in_arb_request),
         .m_arb_enable(p_icm_arb_in_arb_enable),
         .m_arb_read(p_icm_arb_in_arb_read),
         .m_arb_write(p_icm_arb_in_arb_write),
         .m_arb_burstcount(p_icm_arb_in_arb_burstcount),
         .m_arb_address(p_icm_arb_in_arb_address),
         .m_arb_writedata(p_icm_arb_in_arb_writedata),
         .m_arb_byteenable(p_icm_arb_in_arb_byteenable),
         .m_arb_stall(p_icm_arb_in_arb_stall),
         .m_wrp_ack(p_icm_arb_in_wrp_ack),
         .m_rrp_datavalid(p_icm_arb_in_rrp_datavalid),
         .m_rrp_data(p_icm_arb_in_rrp_data),
         // ICM mout
         .mout_arb_request(p_icm_arb_out_arb_request[0]),
         .mout_arb_enable(p_icm_arb_out_arb_enable[0]),
         .mout_arb_read(p_icm_arb_out_arb_read[0]),
         .mout_arb_write(p_icm_arb_out_arb_write[0]),
         .mout_arb_burstcount(p_icm_arb_out_arb_burstcount[0]),
         .mout_arb_address(p_icm_arb_out_arb_address[0]),
         .mout_arb_writedata(p_icm_arb_out_arb_writedata[0]),
         .mout_arb_byteenable(p_icm_arb_out_arb_byteenable[0]),
         .mout_arb_id(),
         .mout_arb_stall(p_icm_arb_out_arb_stall[0]),
         .mout_wrp_ack(p_icm_arb_out_wrp_ack[0]),
         .mout_rrp_datavalid(p_icm_arb_out_rrp_datavalid[0]),
         .mout_rrp_data(p_icm_arb_out_rrp_data[0])
      );

      assign counter_increment_Krnl_IntraE = p_icm_arb_out_arb_address[0];
      assign p_icm_arb_out_rrp_data[0] = counter_result_Krnl_IntraE;
      assign counter_enable_Krnl_IntraE = p_icm_arb_out_arb_read[0];
      assign p_icm_arb_out_rrp_datavalid[0] = counter_resultvalid_Krnl_IntraE;
      assign p_icm_arb_out_arb_stall[0] = counter_stall_Krnl_IntraE;
      assign p_icm_arb_out_wrp_ack[0] = 1'b1;
   end
   endgenerate

   assign counter_resetn_Krnl_Store = (resetn & ~(Krnl_Store_start));
   // INST Krnl_Store_printf_counter of acl_printf_counter
   acl_printf_counter Krnl_Store_printf_counter
   (
      .clock(clock),
      .resetn(counter_resetn_Krnl_Store),
      .enable(counter_enable_Krnl_Store),
      .i_counter_reset(counter_reset_Krnl_Store),
      .i_increment(counter_increment_Krnl_Store),
      .i_init(counter_init_Krnl_Store),
      .i_limit(counter_limit_Krnl_Store),
      .o_size(counter_size_Krnl_Store),
      .o_resultvalid(counter_resultvalid_Krnl_Store),
      .o_result(counter_result_Krnl_Store),
      .o_full(counter_full_Krnl_Store),
      .o_stall(counter_stall_Krnl_Store)
   );

   generate
   begin:printf_system_Krnl_Store
      logic p_icm_arb_in_arb_request [4];
      logic p_icm_arb_in_arb_enable [4];
      logic p_icm_arb_in_arb_read [4];
      logic p_icm_arb_in_arb_write [4];
      logic [5:0] p_icm_arb_in_arb_burstcount [4];
      logic [26:0] p_icm_arb_in_arb_address [4];
      logic [255:0] p_icm_arb_in_arb_writedata [4];
      logic [31:0] p_icm_arb_in_arb_byteenable [4];
      logic p_icm_arb_in_arb_stall [4];
      logic p_icm_arb_in_wrp_ack [4];
      logic p_icm_arb_in_rrp_datavalid [4];
      logic [255:0] p_icm_arb_in_rrp_data [4];
      logic p_icm_arb_out_arb_request [1];
      logic p_icm_arb_out_arb_enable [1];
      logic p_icm_arb_out_arb_read [1];
      logic p_icm_arb_out_arb_write [1];
      logic [5:0] p_icm_arb_out_arb_burstcount [1];
      logic [26:0] p_icm_arb_out_arb_address [1];
      logic [255:0] p_icm_arb_out_arb_writedata [1];
      logic [31:0] p_icm_arb_out_arb_byteenable [1];
      logic p_icm_arb_out_arb_stall [1];
      logic p_icm_arb_out_wrp_ack [1];
      logic p_icm_arb_out_rrp_datavalid [1];
      logic [255:0] p_icm_arb_out_rrp_data [1];

      for( __i = 0; __i < 4; __i = __i + 1 )
      begin:p
         // INST p_avm_to_ic of acl_avm_to_ic
         acl_avm_to_ic
         #(
            .DATA_W(256),
            .WRITEDATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(32),
            .BYTEENA_W(32),
            .ADDR_SHIFT(0)
         )
         p_avm_to_ic
         (
            // AVM avm
            .avm_enable(printf_avm_Krnl_Store_enable[__i]),
            .avm_read(printf_avm_Krnl_Store_read[__i]),
            .avm_write(printf_avm_Krnl_Store_write[__i]),
            .avm_burstcount(printf_avm_Krnl_Store_burstcount[__i]),
            .avm_address(printf_avm_Krnl_Store_address[__i]),
            .avm_writedata(printf_avm_Krnl_Store_writedata[__i]),
            .avm_byteenable(printf_avm_Krnl_Store_byteenable[__i]),
            .avm_waitrequest(printf_avm_Krnl_Store_waitrequest[__i]),
            .avm_readdata(printf_avm_Krnl_Store_readdata[__i]),
            .avm_readdatavalid(printf_avm_Krnl_Store_readdatavalid[__i]),
            // ICM ic
            .ic_arb_request(p_icm_arb_in_arb_request[__i]),
            .ic_arb_enable(p_icm_arb_in_arb_enable[__i]),
            .ic_arb_read(p_icm_arb_in_arb_read[__i]),
            .ic_arb_write(p_icm_arb_in_arb_write[__i]),
            .ic_arb_burstcount(p_icm_arb_in_arb_burstcount[__i]),
            .ic_arb_address(p_icm_arb_in_arb_address[__i]),
            .ic_arb_writedata(p_icm_arb_in_arb_writedata[__i]),
            .ic_arb_byteenable(p_icm_arb_in_arb_byteenable[__i]),
            .ic_arb_stall(p_icm_arb_in_arb_stall[__i]),
            .ic_wrp_ack(p_icm_arb_in_wrp_ack[__i]),
            .ic_rrp_datavalid(p_icm_arb_in_rrp_datavalid[__i]),
            .ic_rrp_data(p_icm_arb_in_rrp_data[__i])
         );

      end

      // INST printf_ic_interconnect of Krnl_GA_system_interconnect_7
      Krnl_GA_system_interconnect_7 printf_ic_interconnect
      (
         .clock(clock),
         .resetn(resetn),
         // ICM m
         .m_arb_request(p_icm_arb_in_arb_request),
         .m_arb_enable(p_icm_arb_in_arb_enable),
         .m_arb_read(p_icm_arb_in_arb_read),
         .m_arb_write(p_icm_arb_in_arb_write),
         .m_arb_burstcount(p_icm_arb_in_arb_burstcount),
         .m_arb_address(p_icm_arb_in_arb_address),
         .m_arb_writedata(p_icm_arb_in_arb_writedata),
         .m_arb_byteenable(p_icm_arb_in_arb_byteenable),
         .m_arb_stall(p_icm_arb_in_arb_stall),
         .m_wrp_ack(p_icm_arb_in_wrp_ack),
         .m_rrp_datavalid(p_icm_arb_in_rrp_datavalid),
         .m_rrp_data(p_icm_arb_in_rrp_data),
         // ICM mout
         .mout_arb_request(p_icm_arb_out_arb_request[0]),
         .mout_arb_enable(p_icm_arb_out_arb_enable[0]),
         .mout_arb_read(p_icm_arb_out_arb_read[0]),
         .mout_arb_write(p_icm_arb_out_arb_write[0]),
         .mout_arb_burstcount(p_icm_arb_out_arb_burstcount[0]),
         .mout_arb_address(p_icm_arb_out_arb_address[0]),
         .mout_arb_writedata(p_icm_arb_out_arb_writedata[0]),
         .mout_arb_byteenable(p_icm_arb_out_arb_byteenable[0]),
         .mout_arb_id(),
         .mout_arb_stall(p_icm_arb_out_arb_stall[0]),
         .mout_wrp_ack(p_icm_arb_out_wrp_ack[0]),
         .mout_rrp_datavalid(p_icm_arb_out_rrp_datavalid[0]),
         .mout_rrp_data(p_icm_arb_out_rrp_data[0])
      );

      assign counter_increment_Krnl_Store = p_icm_arb_out_arb_address[0];
      assign p_icm_arb_out_rrp_data[0] = counter_result_Krnl_Store;
      assign counter_enable_Krnl_Store = p_icm_arb_out_arb_read[0];
      assign p_icm_arb_out_rrp_datavalid[0] = counter_resultvalid_Krnl_Store;
      assign p_icm_arb_out_arb_stall[0] = counter_stall_Krnl_Store;
      assign p_icm_arb_out_wrp_ack[0] = 1'b1;
   end
   endgenerate

endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_Conform_top_wrapper_0
/////////////////////////////////////////////////////////////////
module Krnl_Conform_top_wrapper_0
(
   input logic start,
   input logic [223:0] kernel_arguments,
   input logic [31:0] work_dim,
   input logic [31:0] global_offset [2:0],
   output logic kernel_valid_out,
   output logic has_a_write_pending,
   output logic has_a_lsu_active,
   input logic [31:0] global_id [2:0],
   input logic [31:0] local_id [2:0],
   input logic [31:0] group_id [2:0],
   input logic [31:0] global_size [2:0],
   input logic [31:0] local_size [2:0],
   input logic [31:0] num_groups [2:0],
   input logic [31:0] workgroup_size,
   output logic kernel_stall_out,
   input logic kernel_valid_in,
   input logic clock,
   input logic resetn,
   input logic clock2x,
   // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable,
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read,
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write,
   output logic [4:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount,
   output logic [30:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address,
   output logic [511:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata,
   output logic [63:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest,
   input logic [511:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack,
   // AVM avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0
   output logic avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_enable,
   output logic avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_read,
   output logic avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_write,
   output logic [4:0] avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_burstcount,
   output logic [30:0] avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_address,
   output logic [511:0] avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writedata,
   output logic [63:0] avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_byteenable,
   input logic avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_waitrequest,
   input logic [511:0] avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdata,
   input logic avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid,
   input logic avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writeack,
   // AVM avm_local_bb4_st__inst0
   output logic avm_local_bb4_st__inst0_enable,
   output logic avm_local_bb4_st__inst0_read,
   output logic avm_local_bb4_st__inst0_write,
   output logic [4:0] avm_local_bb4_st__inst0_burstcount,
   output logic [30:0] avm_local_bb4_st__inst0_address,
   output logic [511:0] avm_local_bb4_st__inst0_writedata,
   output logic [63:0] avm_local_bb4_st__inst0_byteenable,
   input logic avm_local_bb4_st__inst0_waitrequest,
   input logic [511:0] avm_local_bb4_st__inst0_readdata,
   input logic avm_local_bb4_st__inst0_readdatavalid,
   input logic avm_local_bb4_st__inst0_writeack,
   // AVM avm_local_bb5_ld__inst0
   output logic avm_local_bb5_ld__inst0_enable,
   output logic avm_local_bb5_ld__inst0_read,
   output logic avm_local_bb5_ld__inst0_write,
   output logic [4:0] avm_local_bb5_ld__inst0_burstcount,
   output logic [30:0] avm_local_bb5_ld__inst0_address,
   output logic [511:0] avm_local_bb5_ld__inst0_writedata,
   output logic [63:0] avm_local_bb5_ld__inst0_byteenable,
   input logic avm_local_bb5_ld__inst0_waitrequest,
   input logic [511:0] avm_local_bb5_ld__inst0_readdata,
   input logic avm_local_bb5_ld__inst0_readdatavalid,
   input logic avm_local_bb5_ld__inst0_writeack,
   // AVM avm_local_bb5_ld__u15_inst0
   output logic avm_local_bb5_ld__u15_inst0_enable,
   output logic avm_local_bb5_ld__u15_inst0_read,
   output logic avm_local_bb5_ld__u15_inst0_write,
   output logic [4:0] avm_local_bb5_ld__u15_inst0_burstcount,
   output logic [30:0] avm_local_bb5_ld__u15_inst0_address,
   output logic [511:0] avm_local_bb5_ld__u15_inst0_writedata,
   output logic [63:0] avm_local_bb5_ld__u15_inst0_byteenable,
   input logic avm_local_bb5_ld__u15_inst0_waitrequest,
   input logic [511:0] avm_local_bb5_ld__u15_inst0_readdata,
   input logic avm_local_bb5_ld__u15_inst0_readdatavalid,
   input logic avm_local_bb5_ld__u15_inst0_writeack,
   // AVM avm_local_bb5_ld__u16_inst0
   output logic avm_local_bb5_ld__u16_inst0_enable,
   output logic avm_local_bb5_ld__u16_inst0_read,
   output logic avm_local_bb5_ld__u16_inst0_write,
   output logic [4:0] avm_local_bb5_ld__u16_inst0_burstcount,
   output logic [30:0] avm_local_bb5_ld__u16_inst0_address,
   output logic [511:0] avm_local_bb5_ld__u16_inst0_writedata,
   output logic [63:0] avm_local_bb5_ld__u16_inst0_byteenable,
   input logic avm_local_bb5_ld__u16_inst0_waitrequest,
   input logic [511:0] avm_local_bb5_ld__u16_inst0_readdata,
   input logic avm_local_bb5_ld__u16_inst0_readdatavalid,
   input logic avm_local_bb5_ld__u16_inst0_writeack,
   // AVM avm_local_bb5_ld__u17_inst0
   output logic avm_local_bb5_ld__u17_inst0_enable,
   output logic avm_local_bb5_ld__u17_inst0_read,
   output logic avm_local_bb5_ld__u17_inst0_write,
   output logic [4:0] avm_local_bb5_ld__u17_inst0_burstcount,
   output logic [30:0] avm_local_bb5_ld__u17_inst0_address,
   output logic [511:0] avm_local_bb5_ld__u17_inst0_writedata,
   output logic [63:0] avm_local_bb5_ld__u17_inst0_byteenable,
   input logic avm_local_bb5_ld__u17_inst0_waitrequest,
   input logic [511:0] avm_local_bb5_ld__u17_inst0_readdata,
   input logic avm_local_bb5_ld__u17_inst0_readdatavalid,
   input logic avm_local_bb5_ld__u17_inst0_writeack,
   // AVM avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0
   output logic avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_enable,
   output logic avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_read,
   output logic avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_write,
   output logic [4:0] avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_burstcount,
   output logic [30:0] avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_address,
   output logic [511:0] avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writedata,
   output logic [63:0] avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_byteenable,
   input logic avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_waitrequest,
   input logic [511:0] avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdata,
   input logic avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdatavalid,
   input logic avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writeack,
   // AVM avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0
   output logic avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_enable,
   output logic avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_read,
   output logic avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_write,
   output logic [4:0] avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_burstcount,
   output logic [30:0] avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_address,
   output logic [511:0] avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writedata,
   output logic [63:0] avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_byteenable,
   input logic avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_waitrequest,
   input logic [511:0] avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdata,
   input logic avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdatavalid,
   input logic avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writeack,
   // AVST avm_channel_id_chan_GA2Conf_genotype_read
   input logic avm_channel_id_chan_GA2Conf_genotype_read_valid,
   output logic avm_channel_id_chan_GA2Conf_genotype_read_ready,
   input logic [31:0] avm_channel_id_chan_GA2Conf_genotype_read_data,
   // AVST avm_channel_id_chan_GA2Conf_active_read
   input logic avm_channel_id_chan_GA2Conf_active_read_valid,
   output logic avm_channel_id_chan_GA2Conf_active_read_ready,
   input logic [7:0] avm_channel_id_chan_GA2Conf_active_read_data,
   // AVST avm_channel_id_chan_GA2Conf_cnt_read
   input logic avm_channel_id_chan_GA2Conf_cnt_read_valid,
   output logic avm_channel_id_chan_GA2Conf_cnt_read_ready,
   input logic [31:0] avm_channel_id_chan_GA2Conf_cnt_read_data,
   // AVST avm_channel_id_chan_GA2Conf_mode_read
   input logic avm_channel_id_chan_GA2Conf_mode_read_valid,
   output logic avm_channel_id_chan_GA2Conf_mode_read_ready,
   input logic [7:0] avm_channel_id_chan_GA2Conf_mode_read_data,
   // AVST avm_channel_id_chan_Conf2Intere_active_write
   output logic avm_channel_id_chan_Conf2Intere_active_write_valid,
   input logic avm_channel_id_chan_Conf2Intere_active_write_ready,
   output logic [7:0] avm_channel_id_chan_Conf2Intere_active_write_data,
   input logic avm_channel_id_chan_Conf2Intere_active_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intere_cnt_write
   output logic avm_channel_id_chan_Conf2Intere_cnt_write_valid,
   input logic avm_channel_id_chan_Conf2Intere_cnt_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intere_cnt_write_data,
   input logic avm_channel_id_chan_Conf2Intere_cnt_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intere_mode_write
   output logic avm_channel_id_chan_Conf2Intere_mode_write_valid,
   input logic avm_channel_id_chan_Conf2Intere_mode_write_ready,
   output logic [7:0] avm_channel_id_chan_Conf2Intere_mode_write_data,
   input logic avm_channel_id_chan_Conf2Intere_mode_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intere_x_write
   output logic avm_channel_id_chan_Conf2Intere_x_write_valid,
   input logic avm_channel_id_chan_Conf2Intere_x_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intere_x_write_data,
   input logic avm_channel_id_chan_Conf2Intere_x_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intere_y_write
   output logic avm_channel_id_chan_Conf2Intere_y_write_valid,
   input logic avm_channel_id_chan_Conf2Intere_y_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intere_y_write_data,
   input logic avm_channel_id_chan_Conf2Intere_y_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intere_z_write
   output logic avm_channel_id_chan_Conf2Intere_z_write_valid,
   input logic avm_channel_id_chan_Conf2Intere_z_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intere_z_write_data,
   input logic avm_channel_id_chan_Conf2Intere_z_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intrae_active_write
   output logic avm_channel_id_chan_Conf2Intrae_active_write_valid,
   input logic avm_channel_id_chan_Conf2Intrae_active_write_ready,
   output logic [7:0] avm_channel_id_chan_Conf2Intrae_active_write_data,
   input logic avm_channel_id_chan_Conf2Intrae_active_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intrae_cnt_write
   output logic avm_channel_id_chan_Conf2Intrae_cnt_write_valid,
   input logic avm_channel_id_chan_Conf2Intrae_cnt_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intrae_cnt_write_data,
   input logic avm_channel_id_chan_Conf2Intrae_cnt_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intrae_mode_write
   output logic avm_channel_id_chan_Conf2Intrae_mode_write_valid,
   input logic avm_channel_id_chan_Conf2Intrae_mode_write_ready,
   output logic [7:0] avm_channel_id_chan_Conf2Intrae_mode_write_data,
   input logic avm_channel_id_chan_Conf2Intrae_mode_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intrae_x_write
   output logic avm_channel_id_chan_Conf2Intrae_x_write_valid,
   input logic avm_channel_id_chan_Conf2Intrae_x_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intrae_x_write_data,
   input logic avm_channel_id_chan_Conf2Intrae_x_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intrae_y_write
   output logic avm_channel_id_chan_Conf2Intrae_y_write_valid,
   input logic avm_channel_id_chan_Conf2Intrae_y_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intrae_y_write_data,
   input logic avm_channel_id_chan_Conf2Intrae_y_write_almostfull,
   // AVST avm_channel_id_chan_Conf2Intrae_z_write
   output logic avm_channel_id_chan_Conf2Intrae_z_write_valid,
   input logic avm_channel_id_chan_Conf2Intrae_z_write_ready,
   output logic [31:0] avm_channel_id_chan_Conf2Intrae_z_write_data,
   input logic avm_channel_id_chan_Conf2Intrae_z_write_almostfull,
   // AVM p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid
);
   genvar __i;
   genvar __j;
   genvar __k;
   genvar __l;
   genvar __m;
   logic lmem_invalid_single_bit;
   logic [3:0] lmem_invalid_aspaces;
   logic local_avm_aspace6_enable [1][6];
   logic local_avm_aspace6_read [1][6];
   logic local_avm_aspace6_write [1][6];
   logic local_avm_aspace6_burstcount [1][6];
   logic [31:0] local_avm_aspace6_address [1][6];
   logic [63:0] local_avm_aspace6_writedata [1][6];
   logic [7:0] local_avm_aspace6_byteenable [1][6];
   logic local_avm_aspace6_waitrequest [1][6];
   logic [63:0] local_avm_aspace6_readdata [1][6];
   logic local_avm_aspace6_readdatavalid [1][6];
   logic local_avm_aspace6_writeack [1][6];
   logic local_avm_aspace7_enable [1][3];
   logic local_avm_aspace7_read [1][3];
   logic local_avm_aspace7_write [1][3];
   logic local_avm_aspace7_burstcount [1][3];
   logic [31:0] local_avm_aspace7_address [1][3];
   logic [31:0] local_avm_aspace7_writedata [1][3];
   logic [3:0] local_avm_aspace7_byteenable [1][3];
   logic local_avm_aspace7_waitrequest [1][3];
   logic [31:0] local_avm_aspace7_readdata [1][3];
   logic local_avm_aspace7_readdatavalid [1][3];
   logic local_avm_aspace7_writeack [1][3];
   logic local_avm_aspace8_enable [1][3];
   logic local_avm_aspace8_read [1][3];
   logic local_avm_aspace8_write [1][3];
   logic local_avm_aspace8_burstcount [1][3];
   logic [31:0] local_avm_aspace8_address [1][3];
   logic [31:0] local_avm_aspace8_writedata [1][3];
   logic [3:0] local_avm_aspace8_byteenable [1][3];
   logic local_avm_aspace8_waitrequest [1][3];
   logic [31:0] local_avm_aspace8_readdata [1][3];
   logic local_avm_aspace8_readdatavalid [1][3];
   logic local_avm_aspace8_writeack [1][3];
   logic local_avm_aspace9_enable [1][3];
   logic local_avm_aspace9_read [1][3];
   logic local_avm_aspace9_write [1][3];
   logic local_avm_aspace9_burstcount [1][3];
   logic [31:0] local_avm_aspace9_address [1][3];
   logic [31:0] local_avm_aspace9_writedata [1][3];
   logic [3:0] local_avm_aspace9_byteenable [1][3];
   logic local_avm_aspace9_waitrequest [1][3];
   logic [31:0] local_avm_aspace9_readdata [1][3];
   logic local_avm_aspace9_readdatavalid [1][3];
   logic local_avm_aspace9_writeack [1][3];

   // INST kernel of Krnl_Conform_function_wrapper
   Krnl_Conform_function_wrapper kernel
   (
      .local_router_hang(lmem_invalid_single_bit),
      .start(start),
      .kernel_arguments(kernel_arguments),
      .work_dim(work_dim),
      .global_offset_0(global_offset[0]),
      .global_offset_1(global_offset[1]),
      .global_offset_2(global_offset[2]),
      .kernel_valid_out(kernel_valid_out),
      .has_a_write_pending(has_a_write_pending),
      .has_a_lsu_active(has_a_lsu_active),
      .global_id_0(global_id[0]),
      .global_id_1(global_id[1]),
      .global_id_2(global_id[2]),
      .local_id_0(local_id[0]),
      .local_id_1(local_id[1]),
      .local_id_2(local_id[2]),
      .group_id_0(group_id[0]),
      .group_id_1(group_id[1]),
      .group_id_2(group_id[2]),
      .global_size_0(global_size[0]),
      .global_size_1(global_size[1]),
      .global_size_2(global_size[2]),
      .local_size_0(local_size[0]),
      .local_size_1(local_size[1]),
      .local_size_2(local_size[2]),
      .num_groups_0(num_groups[0]),
      .num_groups_1(num_groups[1]),
      .num_groups_2(num_groups[2]),
      .workgroup_size(workgroup_size),
      .kernel_stall_out(kernel_stall_out),
      .kernel_valid_in(kernel_valid_in),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack),
      // AVM avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_enable(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_enable),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_read(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_read),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_write(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_write),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_burstcount(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_burstcount),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_address(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_address),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writedata(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writedata),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_byteenable(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_byteenable),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_waitrequest(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_waitrequest),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdata(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdata),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid),
      .avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writeack(avm_local_bb0_ld_memcoalesce_KerConst_load_0_inst0_writeack),
      // AVM avm_local_bb4_st__inst0
      .avm_local_bb4_st__inst0_enable(avm_local_bb4_st__inst0_enable),
      .avm_local_bb4_st__inst0_read(avm_local_bb4_st__inst0_read),
      .avm_local_bb4_st__inst0_write(avm_local_bb4_st__inst0_write),
      .avm_local_bb4_st__inst0_burstcount(avm_local_bb4_st__inst0_burstcount),
      .avm_local_bb4_st__inst0_address(avm_local_bb4_st__inst0_address),
      .avm_local_bb4_st__inst0_writedata(avm_local_bb4_st__inst0_writedata),
      .avm_local_bb4_st__inst0_byteenable(avm_local_bb4_st__inst0_byteenable),
      .avm_local_bb4_st__inst0_waitrequest(avm_local_bb4_st__inst0_waitrequest),
      .avm_local_bb4_st__inst0_readdata(avm_local_bb4_st__inst0_readdata),
      .avm_local_bb4_st__inst0_readdatavalid(avm_local_bb4_st__inst0_readdatavalid),
      .avm_local_bb4_st__inst0_writeack(avm_local_bb4_st__inst0_writeack),
      // AVM avm_local_bb5_ld__inst0
      .avm_local_bb5_ld__inst0_enable(avm_local_bb5_ld__inst0_enable),
      .avm_local_bb5_ld__inst0_read(avm_local_bb5_ld__inst0_read),
      .avm_local_bb5_ld__inst0_write(avm_local_bb5_ld__inst0_write),
      .avm_local_bb5_ld__inst0_burstcount(avm_local_bb5_ld__inst0_burstcount),
      .avm_local_bb5_ld__inst0_address(avm_local_bb5_ld__inst0_address),
      .avm_local_bb5_ld__inst0_writedata(avm_local_bb5_ld__inst0_writedata),
      .avm_local_bb5_ld__inst0_byteenable(avm_local_bb5_ld__inst0_byteenable),
      .avm_local_bb5_ld__inst0_waitrequest(avm_local_bb5_ld__inst0_waitrequest),
      .avm_local_bb5_ld__inst0_readdata(avm_local_bb5_ld__inst0_readdata),
      .avm_local_bb5_ld__inst0_readdatavalid(avm_local_bb5_ld__inst0_readdatavalid),
      .avm_local_bb5_ld__inst0_writeack(avm_local_bb5_ld__inst0_writeack),
      // AVM avm_local_bb5_ld__u15_inst0
      .avm_local_bb5_ld__u15_inst0_enable(avm_local_bb5_ld__u15_inst0_enable),
      .avm_local_bb5_ld__u15_inst0_read(avm_local_bb5_ld__u15_inst0_read),
      .avm_local_bb5_ld__u15_inst0_write(avm_local_bb5_ld__u15_inst0_write),
      .avm_local_bb5_ld__u15_inst0_burstcount(avm_local_bb5_ld__u15_inst0_burstcount),
      .avm_local_bb5_ld__u15_inst0_address(avm_local_bb5_ld__u15_inst0_address),
      .avm_local_bb5_ld__u15_inst0_writedata(avm_local_bb5_ld__u15_inst0_writedata),
      .avm_local_bb5_ld__u15_inst0_byteenable(avm_local_bb5_ld__u15_inst0_byteenable),
      .avm_local_bb5_ld__u15_inst0_waitrequest(avm_local_bb5_ld__u15_inst0_waitrequest),
      .avm_local_bb5_ld__u15_inst0_readdata(avm_local_bb5_ld__u15_inst0_readdata),
      .avm_local_bb5_ld__u15_inst0_readdatavalid(avm_local_bb5_ld__u15_inst0_readdatavalid),
      .avm_local_bb5_ld__u15_inst0_writeack(avm_local_bb5_ld__u15_inst0_writeack),
      // AVM avm_local_bb5_ld__u16_inst0
      .avm_local_bb5_ld__u16_inst0_enable(avm_local_bb5_ld__u16_inst0_enable),
      .avm_local_bb5_ld__u16_inst0_read(avm_local_bb5_ld__u16_inst0_read),
      .avm_local_bb5_ld__u16_inst0_write(avm_local_bb5_ld__u16_inst0_write),
      .avm_local_bb5_ld__u16_inst0_burstcount(avm_local_bb5_ld__u16_inst0_burstcount),
      .avm_local_bb5_ld__u16_inst0_address(avm_local_bb5_ld__u16_inst0_address),
      .avm_local_bb5_ld__u16_inst0_writedata(avm_local_bb5_ld__u16_inst0_writedata),
      .avm_local_bb5_ld__u16_inst0_byteenable(avm_local_bb5_ld__u16_inst0_byteenable),
      .avm_local_bb5_ld__u16_inst0_waitrequest(avm_local_bb5_ld__u16_inst0_waitrequest),
      .avm_local_bb5_ld__u16_inst0_readdata(avm_local_bb5_ld__u16_inst0_readdata),
      .avm_local_bb5_ld__u16_inst0_readdatavalid(avm_local_bb5_ld__u16_inst0_readdatavalid),
      .avm_local_bb5_ld__u16_inst0_writeack(avm_local_bb5_ld__u16_inst0_writeack),
      // AVM avm_local_bb5_ld__u17_inst0
      .avm_local_bb5_ld__u17_inst0_enable(avm_local_bb5_ld__u17_inst0_enable),
      .avm_local_bb5_ld__u17_inst0_read(avm_local_bb5_ld__u17_inst0_read),
      .avm_local_bb5_ld__u17_inst0_write(avm_local_bb5_ld__u17_inst0_write),
      .avm_local_bb5_ld__u17_inst0_burstcount(avm_local_bb5_ld__u17_inst0_burstcount),
      .avm_local_bb5_ld__u17_inst0_address(avm_local_bb5_ld__u17_inst0_address),
      .avm_local_bb5_ld__u17_inst0_writedata(avm_local_bb5_ld__u17_inst0_writedata),
      .avm_local_bb5_ld__u17_inst0_byteenable(avm_local_bb5_ld__u17_inst0_byteenable),
      .avm_local_bb5_ld__u17_inst0_waitrequest(avm_local_bb5_ld__u17_inst0_waitrequest),
      .avm_local_bb5_ld__u17_inst0_readdata(avm_local_bb5_ld__u17_inst0_readdata),
      .avm_local_bb5_ld__u17_inst0_readdatavalid(avm_local_bb5_ld__u17_inst0_readdatavalid),
      .avm_local_bb5_ld__u17_inst0_writeack(avm_local_bb5_ld__u17_inst0_writeack),
      // AVM avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_enable(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_enable),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_read(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_read),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_write(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_write),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_burstcount(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_burstcount),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_address(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_address),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writedata(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writedata),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_byteenable(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_byteenable),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_waitrequest(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_waitrequest),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdata(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdata),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdatavalid(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_readdatavalid),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writeack(avm_local_bb5_ld_memcoalesce_KerConst_load_012_inst0_writeack),
      // AVM avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_enable(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_enable),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_read(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_read),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_write(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_write),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_burstcount(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_burstcount),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_address(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_address),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writedata(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writedata),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_byteenable(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_byteenable),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_waitrequest(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_waitrequest),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdata(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdata),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdatavalid(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_readdatavalid),
      .avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writeack(avm_local_bb5_ld_memcoalesce_KerConst_load_09_inst0_writeack),
      // AVM avm_local_bb3_st__inst0
      .avm_local_bb3_st__inst0_enable(local_avm_aspace6_enable[0][0]),
      .avm_local_bb3_st__inst0_read(local_avm_aspace6_read[0][0]),
      .avm_local_bb3_st__inst0_write(local_avm_aspace6_write[0][0]),
      .avm_local_bb3_st__inst0_burstcount(local_avm_aspace6_burstcount[0][0]),
      .avm_local_bb3_st__inst0_address(local_avm_aspace6_address[0][0]),
      .avm_local_bb3_st__inst0_writedata(local_avm_aspace6_writedata[0][0]),
      .avm_local_bb3_st__inst0_byteenable(local_avm_aspace6_byteenable[0][0]),
      .avm_local_bb3_st__inst0_waitrequest(local_avm_aspace6_waitrequest[0][0]),
      .avm_local_bb3_st__inst0_readdata(local_avm_aspace6_readdata[0][0]),
      .avm_local_bb3_st__inst0_readdatavalid(local_avm_aspace6_readdatavalid[0][0]),
      .avm_local_bb3_st__inst0_writeack(local_avm_aspace6_writeack[0][0]),
      // AVM avm_local_bb4_ld__inst0
      .avm_local_bb4_ld__inst0_enable(local_avm_aspace6_enable[0][1]),
      .avm_local_bb4_ld__inst0_read(local_avm_aspace6_read[0][1]),
      .avm_local_bb4_ld__inst0_write(local_avm_aspace6_write[0][1]),
      .avm_local_bb4_ld__inst0_burstcount(local_avm_aspace6_burstcount[0][1]),
      .avm_local_bb4_ld__inst0_address(local_avm_aspace6_address[0][1]),
      .avm_local_bb4_ld__inst0_writedata(local_avm_aspace6_writedata[0][1]),
      .avm_local_bb4_ld__inst0_byteenable(local_avm_aspace6_byteenable[0][1]),
      .avm_local_bb4_ld__inst0_waitrequest(local_avm_aspace6_waitrequest[0][1]),
      .avm_local_bb4_ld__inst0_readdata(local_avm_aspace6_readdata[0][1]),
      .avm_local_bb4_ld__inst0_readdatavalid(local_avm_aspace6_readdatavalid[0][1]),
      .avm_local_bb4_ld__inst0_writeack(local_avm_aspace6_writeack[0][1]),
      // AVM avm_local_bb4_ld__u6_inst0
      .avm_local_bb4_ld__u6_inst0_enable(local_avm_aspace6_enable[0][2]),
      .avm_local_bb4_ld__u6_inst0_read(local_avm_aspace6_read[0][2]),
      .avm_local_bb4_ld__u6_inst0_write(local_avm_aspace6_write[0][2]),
      .avm_local_bb4_ld__u6_inst0_burstcount(local_avm_aspace6_burstcount[0][2]),
      .avm_local_bb4_ld__u6_inst0_address(local_avm_aspace6_address[0][2]),
      .avm_local_bb4_ld__u6_inst0_writedata(local_avm_aspace6_writedata[0][2]),
      .avm_local_bb4_ld__u6_inst0_byteenable(local_avm_aspace6_byteenable[0][2]),
      .avm_local_bb4_ld__u6_inst0_waitrequest(local_avm_aspace6_waitrequest[0][2]),
      .avm_local_bb4_ld__u6_inst0_readdata(local_avm_aspace6_readdata[0][2]),
      .avm_local_bb4_ld__u6_inst0_readdatavalid(local_avm_aspace6_readdatavalid[0][2]),
      .avm_local_bb4_ld__u6_inst0_writeack(local_avm_aspace6_writeack[0][2]),
      // AVM avm_local_bb4_ld_memcoalesce_null_load_05_inst0
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_enable(local_avm_aspace6_enable[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_read(local_avm_aspace6_read[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_write(local_avm_aspace6_write[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_burstcount(local_avm_aspace6_burstcount[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_address(local_avm_aspace6_address[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_writedata(local_avm_aspace6_writedata[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_byteenable(local_avm_aspace6_byteenable[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_waitrequest(local_avm_aspace6_waitrequest[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_readdata(local_avm_aspace6_readdata[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_readdatavalid(local_avm_aspace6_readdatavalid[0][3]),
      .avm_local_bb4_ld_memcoalesce_null_load_05_inst0_writeack(local_avm_aspace6_writeack[0][3]),
      // AVM avm_local_bb4_ld_memcoalesce_null_load_0_inst0
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_enable(local_avm_aspace6_enable[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_read(local_avm_aspace6_read[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_write(local_avm_aspace6_write[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_burstcount(local_avm_aspace6_burstcount[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_address(local_avm_aspace6_address[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_writedata(local_avm_aspace6_writedata[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_byteenable(local_avm_aspace6_byteenable[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_waitrequest(local_avm_aspace6_waitrequest[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_readdata(local_avm_aspace6_readdata[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_readdatavalid(local_avm_aspace6_readdatavalid[0][4]),
      .avm_local_bb4_ld_memcoalesce_null_load_0_inst0_writeack(local_avm_aspace6_writeack[0][4]),
      // AVM avm_local_bb5_ld__u18_inst0
      .avm_local_bb5_ld__u18_inst0_enable(local_avm_aspace6_enable[0][5]),
      .avm_local_bb5_ld__u18_inst0_read(local_avm_aspace6_read[0][5]),
      .avm_local_bb5_ld__u18_inst0_write(local_avm_aspace6_write[0][5]),
      .avm_local_bb5_ld__u18_inst0_burstcount(local_avm_aspace6_burstcount[0][5]),
      .avm_local_bb5_ld__u18_inst0_address(local_avm_aspace6_address[0][5]),
      .avm_local_bb5_ld__u18_inst0_writedata(local_avm_aspace6_writedata[0][5]),
      .avm_local_bb5_ld__u18_inst0_byteenable(local_avm_aspace6_byteenable[0][5]),
      .avm_local_bb5_ld__u18_inst0_waitrequest(local_avm_aspace6_waitrequest[0][5]),
      .avm_local_bb5_ld__u18_inst0_readdata(local_avm_aspace6_readdata[0][5]),
      .avm_local_bb5_ld__u18_inst0_readdatavalid(local_avm_aspace6_readdatavalid[0][5]),
      .avm_local_bb5_ld__u18_inst0_writeack(local_avm_aspace6_writeack[0][5]),
      // AVM avm_local_bb5_ld__u21_inst0
      .avm_local_bb5_ld__u21_inst0_enable(local_avm_aspace7_enable[0][0]),
      .avm_local_bb5_ld__u21_inst0_read(local_avm_aspace7_read[0][0]),
      .avm_local_bb5_ld__u21_inst0_write(local_avm_aspace7_write[0][0]),
      .avm_local_bb5_ld__u21_inst0_burstcount(local_avm_aspace7_burstcount[0][0]),
      .avm_local_bb5_ld__u21_inst0_address(local_avm_aspace7_address[0][0]),
      .avm_local_bb5_ld__u21_inst0_writedata(local_avm_aspace7_writedata[0][0]),
      .avm_local_bb5_ld__u21_inst0_byteenable(local_avm_aspace7_byteenable[0][0]),
      .avm_local_bb5_ld__u21_inst0_waitrequest(local_avm_aspace7_waitrequest[0][0]),
      .avm_local_bb5_ld__u21_inst0_readdata(local_avm_aspace7_readdata[0][0]),
      .avm_local_bb5_ld__u21_inst0_readdatavalid(local_avm_aspace7_readdatavalid[0][0]),
      .avm_local_bb5_ld__u21_inst0_writeack(local_avm_aspace7_writeack[0][0]),
      // AVM avm_local_bb5_st_add241_hfp_inst0
      .avm_local_bb5_st_add241_hfp_inst0_enable(local_avm_aspace7_enable[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_read(local_avm_aspace7_read[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_write(local_avm_aspace7_write[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_burstcount(local_avm_aspace7_burstcount[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_address(local_avm_aspace7_address[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_writedata(local_avm_aspace7_writedata[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_byteenable(local_avm_aspace7_byteenable[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_waitrequest(local_avm_aspace7_waitrequest[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_readdata(local_avm_aspace7_readdata[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_readdatavalid(local_avm_aspace7_readdatavalid[0][1]),
      .avm_local_bb5_st_add241_hfp_inst0_writeack(local_avm_aspace7_writeack[0][1]),
      // AVM avm_local_bb6_ld__inst0
      .avm_local_bb6_ld__inst0_enable(local_avm_aspace7_enable[0][2]),
      .avm_local_bb6_ld__inst0_read(local_avm_aspace7_read[0][2]),
      .avm_local_bb6_ld__inst0_write(local_avm_aspace7_write[0][2]),
      .avm_local_bb6_ld__inst0_burstcount(local_avm_aspace7_burstcount[0][2]),
      .avm_local_bb6_ld__inst0_address(local_avm_aspace7_address[0][2]),
      .avm_local_bb6_ld__inst0_writedata(local_avm_aspace7_writedata[0][2]),
      .avm_local_bb6_ld__inst0_byteenable(local_avm_aspace7_byteenable[0][2]),
      .avm_local_bb6_ld__inst0_waitrequest(local_avm_aspace7_waitrequest[0][2]),
      .avm_local_bb6_ld__inst0_readdata(local_avm_aspace7_readdata[0][2]),
      .avm_local_bb6_ld__inst0_readdatavalid(local_avm_aspace7_readdatavalid[0][2]),
      .avm_local_bb6_ld__inst0_writeack(local_avm_aspace7_writeack[0][2]),
      // AVM avm_local_bb5_ld__u22_inst0
      .avm_local_bb5_ld__u22_inst0_enable(local_avm_aspace8_enable[0][0]),
      .avm_local_bb5_ld__u22_inst0_read(local_avm_aspace8_read[0][0]),
      .avm_local_bb5_ld__u22_inst0_write(local_avm_aspace8_write[0][0]),
      .avm_local_bb5_ld__u22_inst0_burstcount(local_avm_aspace8_burstcount[0][0]),
      .avm_local_bb5_ld__u22_inst0_address(local_avm_aspace8_address[0][0]),
      .avm_local_bb5_ld__u22_inst0_writedata(local_avm_aspace8_writedata[0][0]),
      .avm_local_bb5_ld__u22_inst0_byteenable(local_avm_aspace8_byteenable[0][0]),
      .avm_local_bb5_ld__u22_inst0_waitrequest(local_avm_aspace8_waitrequest[0][0]),
      .avm_local_bb5_ld__u22_inst0_readdata(local_avm_aspace8_readdata[0][0]),
      .avm_local_bb5_ld__u22_inst0_readdatavalid(local_avm_aspace8_readdatavalid[0][0]),
      .avm_local_bb5_ld__u22_inst0_writeack(local_avm_aspace8_writeack[0][0]),
      // AVM avm_local_bb5_st_add246_hfp_inst0
      .avm_local_bb5_st_add246_hfp_inst0_enable(local_avm_aspace8_enable[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_read(local_avm_aspace8_read[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_write(local_avm_aspace8_write[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_burstcount(local_avm_aspace8_burstcount[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_address(local_avm_aspace8_address[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_writedata(local_avm_aspace8_writedata[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_byteenable(local_avm_aspace8_byteenable[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_waitrequest(local_avm_aspace8_waitrequest[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_readdata(local_avm_aspace8_readdata[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_readdatavalid(local_avm_aspace8_readdatavalid[0][1]),
      .avm_local_bb5_st_add246_hfp_inst0_writeack(local_avm_aspace8_writeack[0][1]),
      // AVM avm_local_bb6_ld__u25_inst0
      .avm_local_bb6_ld__u25_inst0_enable(local_avm_aspace8_enable[0][2]),
      .avm_local_bb6_ld__u25_inst0_read(local_avm_aspace8_read[0][2]),
      .avm_local_bb6_ld__u25_inst0_write(local_avm_aspace8_write[0][2]),
      .avm_local_bb6_ld__u25_inst0_burstcount(local_avm_aspace8_burstcount[0][2]),
      .avm_local_bb6_ld__u25_inst0_address(local_avm_aspace8_address[0][2]),
      .avm_local_bb6_ld__u25_inst0_writedata(local_avm_aspace8_writedata[0][2]),
      .avm_local_bb6_ld__u25_inst0_byteenable(local_avm_aspace8_byteenable[0][2]),
      .avm_local_bb6_ld__u25_inst0_waitrequest(local_avm_aspace8_waitrequest[0][2]),
      .avm_local_bb6_ld__u25_inst0_readdata(local_avm_aspace8_readdata[0][2]),
      .avm_local_bb6_ld__u25_inst0_readdatavalid(local_avm_aspace8_readdatavalid[0][2]),
      .avm_local_bb6_ld__u25_inst0_writeack(local_avm_aspace8_writeack[0][2]),
      // AVM avm_local_bb5_ld__u23_inst0
      .avm_local_bb5_ld__u23_inst0_enable(local_avm_aspace9_enable[0][0]),
      .avm_local_bb5_ld__u23_inst0_read(local_avm_aspace9_read[0][0]),
      .avm_local_bb5_ld__u23_inst0_write(local_avm_aspace9_write[0][0]),
      .avm_local_bb5_ld__u23_inst0_burstcount(local_avm_aspace9_burstcount[0][0]),
      .avm_local_bb5_ld__u23_inst0_address(local_avm_aspace9_address[0][0]),
      .avm_local_bb5_ld__u23_inst0_writedata(local_avm_aspace9_writedata[0][0]),
      .avm_local_bb5_ld__u23_inst0_byteenable(local_avm_aspace9_byteenable[0][0]),
      .avm_local_bb5_ld__u23_inst0_waitrequest(local_avm_aspace9_waitrequest[0][0]),
      .avm_local_bb5_ld__u23_inst0_readdata(local_avm_aspace9_readdata[0][0]),
      .avm_local_bb5_ld__u23_inst0_readdatavalid(local_avm_aspace9_readdatavalid[0][0]),
      .avm_local_bb5_ld__u23_inst0_writeack(local_avm_aspace9_writeack[0][0]),
      // AVM avm_local_bb5_st_add251_hfp_inst0
      .avm_local_bb5_st_add251_hfp_inst0_enable(local_avm_aspace9_enable[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_read(local_avm_aspace9_read[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_write(local_avm_aspace9_write[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_burstcount(local_avm_aspace9_burstcount[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_address(local_avm_aspace9_address[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_writedata(local_avm_aspace9_writedata[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_byteenable(local_avm_aspace9_byteenable[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_waitrequest(local_avm_aspace9_waitrequest[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_readdata(local_avm_aspace9_readdata[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_readdatavalid(local_avm_aspace9_readdatavalid[0][1]),
      .avm_local_bb5_st_add251_hfp_inst0_writeack(local_avm_aspace9_writeack[0][1]),
      // AVM avm_local_bb6_ld__u27_inst0
      .avm_local_bb6_ld__u27_inst0_enable(local_avm_aspace9_enable[0][2]),
      .avm_local_bb6_ld__u27_inst0_read(local_avm_aspace9_read[0][2]),
      .avm_local_bb6_ld__u27_inst0_write(local_avm_aspace9_write[0][2]),
      .avm_local_bb6_ld__u27_inst0_burstcount(local_avm_aspace9_burstcount[0][2]),
      .avm_local_bb6_ld__u27_inst0_address(local_avm_aspace9_address[0][2]),
      .avm_local_bb6_ld__u27_inst0_writedata(local_avm_aspace9_writedata[0][2]),
      .avm_local_bb6_ld__u27_inst0_byteenable(local_avm_aspace9_byteenable[0][2]),
      .avm_local_bb6_ld__u27_inst0_waitrequest(local_avm_aspace9_waitrequest[0][2]),
      .avm_local_bb6_ld__u27_inst0_readdata(local_avm_aspace9_readdata[0][2]),
      .avm_local_bb6_ld__u27_inst0_readdatavalid(local_avm_aspace9_readdatavalid[0][2]),
      .avm_local_bb6_ld__u27_inst0_writeack(local_avm_aspace9_writeack[0][2]),
      // AVST avst_local_bb3__chan_GA2Conf_genotype_inst0
      .avst_local_bb3__chan_GA2Conf_genotype_inst0_valid(avm_channel_id_chan_GA2Conf_genotype_read_valid),
      .avst_local_bb3__chan_GA2Conf_genotype_inst0_ready(avm_channel_id_chan_GA2Conf_genotype_read_ready),
      .avst_local_bb3__chan_GA2Conf_genotype_inst0_data(avm_channel_id_chan_GA2Conf_genotype_read_data),
      // AVST avst_local_bb4__chan_GA2Conf_active_inst0
      .avst_local_bb4__chan_GA2Conf_active_inst0_valid(avm_channel_id_chan_GA2Conf_active_read_valid),
      .avst_local_bb4__chan_GA2Conf_active_inst0_ready(avm_channel_id_chan_GA2Conf_active_read_ready),
      .avst_local_bb4__chan_GA2Conf_active_inst0_data(avm_channel_id_chan_GA2Conf_active_read_data),
      // AVST avst_local_bb4__chan_GA2Conf_cnt_inst0
      .avst_local_bb4__chan_GA2Conf_cnt_inst0_valid(avm_channel_id_chan_GA2Conf_cnt_read_valid),
      .avst_local_bb4__chan_GA2Conf_cnt_inst0_ready(avm_channel_id_chan_GA2Conf_cnt_read_ready),
      .avst_local_bb4__chan_GA2Conf_cnt_inst0_data(avm_channel_id_chan_GA2Conf_cnt_read_data),
      // AVST avst_local_bb4__chan_GA2Conf_mode_inst0
      .avst_local_bb4__chan_GA2Conf_mode_inst0_valid(avm_channel_id_chan_GA2Conf_mode_read_valid),
      .avst_local_bb4__chan_GA2Conf_mode_inst0_ready(avm_channel_id_chan_GA2Conf_mode_read_ready),
      .avst_local_bb4__chan_GA2Conf_mode_inst0_data(avm_channel_id_chan_GA2Conf_mode_read_data),
      // AVST avst_local_bb6__chan_Conf2Intere_active_inst0
      .avst_local_bb6__chan_Conf2Intere_active_inst0_valid(avm_channel_id_chan_Conf2Intere_active_write_valid),
      .avst_local_bb6__chan_Conf2Intere_active_inst0_ready(avm_channel_id_chan_Conf2Intere_active_write_ready),
      .avst_local_bb6__chan_Conf2Intere_active_inst0_data(avm_channel_id_chan_Conf2Intere_active_write_data),
      .avst_local_bb6__chan_Conf2Intere_active_inst0_almostfull(avm_channel_id_chan_Conf2Intere_active_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intere_cnt_inst0
      .avst_local_bb6__chan_Conf2Intere_cnt_inst0_valid(avm_channel_id_chan_Conf2Intere_cnt_write_valid),
      .avst_local_bb6__chan_Conf2Intere_cnt_inst0_ready(avm_channel_id_chan_Conf2Intere_cnt_write_ready),
      .avst_local_bb6__chan_Conf2Intere_cnt_inst0_data(avm_channel_id_chan_Conf2Intere_cnt_write_data),
      .avst_local_bb6__chan_Conf2Intere_cnt_inst0_almostfull(avm_channel_id_chan_Conf2Intere_cnt_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intere_mode_inst0
      .avst_local_bb6__chan_Conf2Intere_mode_inst0_valid(avm_channel_id_chan_Conf2Intere_mode_write_valid),
      .avst_local_bb6__chan_Conf2Intere_mode_inst0_ready(avm_channel_id_chan_Conf2Intere_mode_write_ready),
      .avst_local_bb6__chan_Conf2Intere_mode_inst0_data(avm_channel_id_chan_Conf2Intere_mode_write_data),
      .avst_local_bb6__chan_Conf2Intere_mode_inst0_almostfull(avm_channel_id_chan_Conf2Intere_mode_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intere_x_inst0
      .avst_local_bb6__chan_Conf2Intere_x_inst0_valid(avm_channel_id_chan_Conf2Intere_x_write_valid),
      .avst_local_bb6__chan_Conf2Intere_x_inst0_ready(avm_channel_id_chan_Conf2Intere_x_write_ready),
      .avst_local_bb6__chan_Conf2Intere_x_inst0_data(avm_channel_id_chan_Conf2Intere_x_write_data),
      .avst_local_bb6__chan_Conf2Intere_x_inst0_almostfull(avm_channel_id_chan_Conf2Intere_x_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intere_y_inst0
      .avst_local_bb6__chan_Conf2Intere_y_inst0_valid(avm_channel_id_chan_Conf2Intere_y_write_valid),
      .avst_local_bb6__chan_Conf2Intere_y_inst0_ready(avm_channel_id_chan_Conf2Intere_y_write_ready),
      .avst_local_bb6__chan_Conf2Intere_y_inst0_data(avm_channel_id_chan_Conf2Intere_y_write_data),
      .avst_local_bb6__chan_Conf2Intere_y_inst0_almostfull(avm_channel_id_chan_Conf2Intere_y_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intere_z_inst0
      .avst_local_bb6__chan_Conf2Intere_z_inst0_valid(avm_channel_id_chan_Conf2Intere_z_write_valid),
      .avst_local_bb6__chan_Conf2Intere_z_inst0_ready(avm_channel_id_chan_Conf2Intere_z_write_ready),
      .avst_local_bb6__chan_Conf2Intere_z_inst0_data(avm_channel_id_chan_Conf2Intere_z_write_data),
      .avst_local_bb6__chan_Conf2Intere_z_inst0_almostfull(avm_channel_id_chan_Conf2Intere_z_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intrae_active_inst0
      .avst_local_bb6__chan_Conf2Intrae_active_inst0_valid(avm_channel_id_chan_Conf2Intrae_active_write_valid),
      .avst_local_bb6__chan_Conf2Intrae_active_inst0_ready(avm_channel_id_chan_Conf2Intrae_active_write_ready),
      .avst_local_bb6__chan_Conf2Intrae_active_inst0_data(avm_channel_id_chan_Conf2Intrae_active_write_data),
      .avst_local_bb6__chan_Conf2Intrae_active_inst0_almostfull(avm_channel_id_chan_Conf2Intrae_active_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intrae_cnt_inst0
      .avst_local_bb6__chan_Conf2Intrae_cnt_inst0_valid(avm_channel_id_chan_Conf2Intrae_cnt_write_valid),
      .avst_local_bb6__chan_Conf2Intrae_cnt_inst0_ready(avm_channel_id_chan_Conf2Intrae_cnt_write_ready),
      .avst_local_bb6__chan_Conf2Intrae_cnt_inst0_data(avm_channel_id_chan_Conf2Intrae_cnt_write_data),
      .avst_local_bb6__chan_Conf2Intrae_cnt_inst0_almostfull(avm_channel_id_chan_Conf2Intrae_cnt_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intrae_mode_inst0
      .avst_local_bb6__chan_Conf2Intrae_mode_inst0_valid(avm_channel_id_chan_Conf2Intrae_mode_write_valid),
      .avst_local_bb6__chan_Conf2Intrae_mode_inst0_ready(avm_channel_id_chan_Conf2Intrae_mode_write_ready),
      .avst_local_bb6__chan_Conf2Intrae_mode_inst0_data(avm_channel_id_chan_Conf2Intrae_mode_write_data),
      .avst_local_bb6__chan_Conf2Intrae_mode_inst0_almostfull(avm_channel_id_chan_Conf2Intrae_mode_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intrae_x_inst0
      .avst_local_bb6__chan_Conf2Intrae_x_inst0_valid(avm_channel_id_chan_Conf2Intrae_x_write_valid),
      .avst_local_bb6__chan_Conf2Intrae_x_inst0_ready(avm_channel_id_chan_Conf2Intrae_x_write_ready),
      .avst_local_bb6__chan_Conf2Intrae_x_inst0_data(avm_channel_id_chan_Conf2Intrae_x_write_data),
      .avst_local_bb6__chan_Conf2Intrae_x_inst0_almostfull(avm_channel_id_chan_Conf2Intrae_x_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intrae_y_inst0
      .avst_local_bb6__chan_Conf2Intrae_y_inst0_valid(avm_channel_id_chan_Conf2Intrae_y_write_valid),
      .avst_local_bb6__chan_Conf2Intrae_y_inst0_ready(avm_channel_id_chan_Conf2Intrae_y_write_ready),
      .avst_local_bb6__chan_Conf2Intrae_y_inst0_data(avm_channel_id_chan_Conf2Intrae_y_write_data),
      .avst_local_bb6__chan_Conf2Intrae_y_inst0_almostfull(avm_channel_id_chan_Conf2Intrae_y_write_almostfull),
      // AVST avst_local_bb6__chan_Conf2Intrae_z_inst0
      .avst_local_bb6__chan_Conf2Intrae_z_inst0_valid(avm_channel_id_chan_Conf2Intrae_z_write_valid),
      .avst_local_bb6__chan_Conf2Intrae_z_inst0_ready(avm_channel_id_chan_Conf2Intrae_z_write_ready),
      .avst_local_bb6__chan_Conf2Intrae_z_inst0_data(avm_channel_id_chan_Conf2Intrae_z_write_data),
      .avst_local_bb6__chan_Conf2Intrae_z_inst0_almostfull(avm_channel_id_chan_Conf2Intrae_z_write_almostfull),
      // AVM p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_read(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_write(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_address(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb4_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid)
   );

   assign lmem_invalid_single_bit = |lmem_invalid_aspaces;
   generate
   begin:local_mem_system_aspace6
      logic local_icm_arb_request [1][6];
      logic local_icm_arb_enable [1][6];
      logic local_icm_arb_read [1][6];
      logic local_icm_arb_write [1][6];
      logic local_icm_arb_burstcount [1][6];
      logic [4:0] local_icm_arb_address [1][6];
      logic [63:0] local_icm_arb_writedata [1][6];
      logic [7:0] local_icm_arb_byteenable [1][6];
      logic local_icm_arb_stall [1][6];
      logic local_icm_wrp_ack [1][6];
      logic local_icm_rrp_datavalid [1][6];
      logic [63:0] local_icm_rrp_data [1][6];
      logic invalid_access_grps;

      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:local_mem_group
         logic [3:0] invalid_access_terms;

         for( __j = 0; __j < 6; __j = __j + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(64),
               .WRITEDATA_W(64),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(8)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace6_enable[__i][__j]),
               .avm_read(local_avm_aspace6_read[__i][__j]),
               .avm_write(local_avm_aspace6_write[__i][__j]),
               .avm_burstcount(local_avm_aspace6_burstcount[__i][__j]),
               .avm_address(local_avm_aspace6_address[__i][__j]),
               .avm_writedata(local_avm_aspace6_writedata[__i][__j]),
               .avm_byteenable(local_avm_aspace6_byteenable[__i][__j]),
               .avm_waitrequest(local_avm_aspace6_waitrequest[__i][__j]),
               .avm_readdata(local_avm_aspace6_readdata[__i][__j]),
               .avm_readdatavalid(local_avm_aspace6_readdatavalid[__i][__j]),
               .avm_writeack(local_avm_aspace6_writeack[__i][__j]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__i][__j]),
               .ic_arb_enable(local_icm_arb_enable[__i][__j]),
               .ic_arb_read(local_icm_arb_read[__i][__j]),
               .ic_arb_write(local_icm_arb_write[__i][__j]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__i][__j]),
               .ic_arb_address(local_icm_arb_address[__i][__j]),
               .ic_arb_writedata(local_icm_arb_writedata[__i][__j]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__i][__j]),
               .ic_arb_stall(local_icm_arb_stall[__i][__j]),
               .ic_wrp_ack(local_icm_wrp_ack[__i][__j]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__i][__j]),
               .ic_rrp_data(local_icm_rrp_data[__i][__j])
            );

         end

         for( __j = 0; __j < 2; __j = __j + 1 )
         begin:bank
            logic port_enable [1:4];
            logic port_read [1:4];
            logic port_write [1:4];
            logic [3:0] port_address [1:4];
            logic [63:0] port_writedata [1:4];
            logic [7:0] port_byteenable [1:4];
            logic port_waitrequest [1:4];
            logic [63:0] port_readdata [1:4];
            logic port_readdatavalid [1:4];

            // INST mem0 of acl_mem2x
            acl_mem2x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(16),
               .WIDTH(64),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("BIDIR_DUAL_PORT"),
               .PREFERRED_WIDTH(160),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .clk2x(clock2x),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2]),
               // AVS avs_port3
               .avs_port3_enable(port_enable[3]),
               .avs_port3_read(port_read[3]),
               .avs_port3_write(port_write[3]),
               .avs_port3_address(port_address[3]),
               .avs_port3_writedata(port_writedata[3]),
               .avs_port3_byteenable(port_byteenable[3]),
               .avs_port3_waitrequest(port_waitrequest[3]),
               .avs_port3_readdata(port_readdata[3]),
               .avs_port3_readdatavalid(port_readdatavalid[3]),
               // AVS avs_port4
               .avs_port4_enable(port_enable[4]),
               .avs_port4_read(port_read[4]),
               .avs_port4_write(port_write[4]),
               .avs_port4_address(port_address[4]),
               .avs_port4_writedata(port_writedata[4]),
               .avs_port4_byteenable(port_byteenable[4]),
               .avs_port4_waitrequest(port_waitrequest[4]),
               .avs_port4_readdata(port_readdata[4]),
               .avs_port4_readdatavalid(port_readdatavalid[4])
            );

         end

         for( __j = 0; __j < 6; __j = __j + 1 )
         begin:router
            logic b_arb_request [2];
            logic b_arb_enable [2];
            logic b_arb_read [2];
            logic b_arb_write [2];
            logic b_arb_burstcount [2];
            logic [3:0] b_arb_address [2];
            logic [63:0] b_arb_writedata [2];
            logic [7:0] b_arb_byteenable [2];
            logic b_arb_stall [2];
            logic b_wrp_ack [2];
            logic b_rrp_datavalid [2];
            logic [63:0] b_rrp_data [2];
            logic [1:0] bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(64),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(5),
               .BYTEENA_W(8),
               .NUM_BANKS(2)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__i][__j]),
               .m_arb_enable(local_icm_arb_enable[__i][__j]),
               .m_arb_read(local_icm_arb_read[__i][__j]),
               .m_arb_write(local_icm_arb_write[__i][__j]),
               .m_arb_burstcount(local_icm_arb_burstcount[__i][__j]),
               .m_arb_address(local_icm_arb_address[__i][__j]),
               .m_arb_writedata(local_icm_arb_writedata[__i][__j]),
               .m_arb_byteenable(local_icm_arb_byteenable[__i][__j]),
               .m_arb_stall(local_icm_arb_stall[__i][__j]),
               .m_wrp_ack(local_icm_wrp_ack[__i][__j]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__i][__j]),
               .m_rrp_data(local_icm_rrp_data[__i][__j]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select[0] = (local_icm_arb_address[__i][__j][4] == 1'b0);
            assign bank_select[1] = (local_icm_arb_address[__i][__j][4] == 1'b1);
         end

         assign invalid_access_grps = |invalid_access_terms;
         // INST acl_ic_local_mem_router_terminator_m1b0 of acl_ic_local_mem_router_terminator
         acl_ic_local_mem_router_terminator
         #(
            .DATA_W(64)
         )
         acl_ic_local_mem_router_terminator_m1b0
         (
            .clock(clock),
            .resetn(resetn),
            .b_arb_request(router[1].b_arb_request[0]),
            .b_arb_read(router[1].b_arb_read[0]),
            .b_arb_write(router[1].b_arb_write[0]),
            .b_arb_stall(router[1].b_arb_stall[0]),
            .b_wrp_ack(router[1].b_wrp_ack[0]),
            .b_rrp_datavalid(router[1].b_rrp_datavalid[0]),
            .b_rrp_data(router[1].b_rrp_data[0]),
            .b_invalid_access(invalid_access_terms[0])
         );

         // INST acl_ic_local_mem_router_terminator_m2b0 of acl_ic_local_mem_router_terminator
         acl_ic_local_mem_router_terminator
         #(
            .DATA_W(64)
         )
         acl_ic_local_mem_router_terminator_m2b0
         (
            .clock(clock),
            .resetn(resetn),
            .b_arb_request(router[2].b_arb_request[0]),
            .b_arb_read(router[2].b_arb_read[0]),
            .b_arb_write(router[2].b_arb_write[0]),
            .b_arb_stall(router[2].b_arb_stall[0]),
            .b_wrp_ack(router[2].b_wrp_ack[0]),
            .b_rrp_datavalid(router[2].b_rrp_datavalid[0]),
            .b_rrp_data(router[2].b_rrp_data[0]),
            .b_invalid_access(invalid_access_terms[1])
         );

         // INST acl_ic_local_mem_router_terminator_m3b1 of acl_ic_local_mem_router_terminator
         acl_ic_local_mem_router_terminator
         #(
            .DATA_W(64)
         )
         acl_ic_local_mem_router_terminator_m3b1
         (
            .clock(clock),
            .resetn(resetn),
            .b_arb_request(router[3].b_arb_request[1]),
            .b_arb_read(router[3].b_arb_read[1]),
            .b_arb_write(router[3].b_arb_write[1]),
            .b_arb_stall(router[3].b_arb_stall[1]),
            .b_wrp_ack(router[3].b_wrp_ack[1]),
            .b_rrp_datavalid(router[3].b_rrp_datavalid[1]),
            .b_rrp_data(router[3].b_rrp_data[1]),
            .b_invalid_access(invalid_access_terms[2])
         );

         // INST acl_ic_local_mem_router_terminator_m4b1 of acl_ic_local_mem_router_terminator
         acl_ic_local_mem_router_terminator
         #(
            .DATA_W(64)
         )
         acl_ic_local_mem_router_terminator_m4b1
         (
            .clock(clock),
            .resetn(resetn),
            .b_arb_request(router[4].b_arb_request[1]),
            .b_arb_read(router[4].b_arb_read[1]),
            .b_arb_write(router[4].b_arb_write[1]),
            .b_arb_stall(router[4].b_arb_stall[1]),
            .b_wrp_ack(router[4].b_wrp_ack[1]),
            .b_rrp_datavalid(router[4].b_rrp_datavalid[1]),
            .b_rrp_data(router[4].b_rrp_data[1]),
            .b_invalid_access(invalid_access_terms[3])
         );

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_0
            Krnl_GA_system_interconnect_0 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port1bank1
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[1];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[1];
            assign icm_in_arb_read[0] = router[0].b_arb_read[1];
            assign icm_in_arb_write[0] = router[0].b_arb_write[1];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[1];
            assign icm_in_arb_address[0] = router[0].b_arb_address[1];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[1];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[1];
            assign router[0].b_arb_stall[1] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[1] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[1] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[1] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_0
            Krnl_GA_system_interconnect_0 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[1].port_enable[1] = icm_out_arb_enable;
            assign bank[1].port_read[1] = icm_out_arb_read;
            assign bank[1].port_write[1] = icm_out_arb_write;
            assign bank[1].port_address[1] = icm_out_arb_address;
            assign bank[1].port_writedata[1] = icm_out_arb_writedata;
            assign bank[1].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[1].port_waitrequest[1];
            assign icm_out_rrp_data = bank[1].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[1].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[4].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[4].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[4].b_arb_read[0];
            assign icm_in_arb_write[0] = router[4].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[4].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[4].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[4].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[4].b_arb_byteenable[0];
            assign router[4].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[4].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[4].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[4].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_1
            Krnl_GA_system_interconnect_1 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port2bank1
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[1];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[1];
            assign icm_in_arb_read[0] = router[1].b_arb_read[1];
            assign icm_in_arb_write[0] = router[1].b_arb_write[1];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[1];
            assign icm_in_arb_address[0] = router[1].b_arb_address[1];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[1];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[1];
            assign router[1].b_arb_stall[1] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[1] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[1] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[1] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_1
            Krnl_GA_system_interconnect_1 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[1].port_enable[2] = icm_out_arb_enable;
            assign bank[1].port_read[2] = icm_out_arb_read;
            assign bank[1].port_write[2] = icm_out_arb_write;
            assign bank[1].port_address[2] = icm_out_arb_address;
            assign bank[1].port_writedata[2] = icm_out_arb_writedata;
            assign bank[1].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[1].port_waitrequest[2];
            assign icm_out_rrp_data = bank[1].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[1].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port3bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[3].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[3].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[3].b_arb_read[0];
            assign icm_in_arb_write[0] = router[3].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[3].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[3].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[3].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[3].b_arb_byteenable[0];
            assign router[3].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[3].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[3].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[3].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_1
            Krnl_GA_system_interconnect_1 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[3] = icm_out_arb_enable;
            assign bank[0].port_read[3] = icm_out_arb_read;
            assign bank[0].port_write[3] = icm_out_arb_write;
            assign bank[0].port_address[3] = icm_out_arb_address;
            assign bank[0].port_writedata[3] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[3];
            assign icm_out_rrp_data = bank[0].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port3bank1
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[2].b_arb_request[1];
            assign icm_in_arb_enable[0] = router[2].b_arb_enable[1];
            assign icm_in_arb_read[0] = router[2].b_arb_read[1];
            assign icm_in_arb_write[0] = router[2].b_arb_write[1];
            assign icm_in_arb_burstcount[0] = router[2].b_arb_burstcount[1];
            assign icm_in_arb_address[0] = router[2].b_arb_address[1];
            assign icm_in_arb_writedata[0] = router[2].b_arb_writedata[1];
            assign icm_in_arb_byteenable[0] = router[2].b_arb_byteenable[1];
            assign router[2].b_arb_stall[1] = icm_in_arb_stall[0];
            assign router[2].b_wrp_ack[1] = icm_in_wrp_ack[0];
            assign router[2].b_rrp_datavalid[1] = icm_in_rrp_datavalid[0];
            assign router[2].b_rrp_data[1] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_1
            Krnl_GA_system_interconnect_1 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[1].port_enable[3] = icm_out_arb_enable;
            assign bank[1].port_read[3] = icm_out_arb_read;
            assign bank[1].port_write[3] = icm_out_arb_write;
            assign bank[1].port_address[3] = icm_out_arb_address;
            assign bank[1].port_writedata[3] = icm_out_arb_writedata;
            assign bank[1].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[1].port_waitrequest[3];
            assign icm_out_rrp_data = bank[1].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[1].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port4bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[5].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[5].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[5].b_arb_read[0];
            assign icm_in_arb_write[0] = router[5].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[5].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[5].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[5].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[5].b_arb_byteenable[0];
            assign router[5].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[5].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[5].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[5].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_1
            Krnl_GA_system_interconnect_1 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[4] = icm_out_arb_enable;
            assign bank[0].port_read[4] = icm_out_arb_read;
            assign bank[0].port_write[4] = icm_out_arb_write;
            assign bank[0].port_address[4] = icm_out_arb_address;
            assign bank[0].port_writedata[4] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[4] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[4];
            assign icm_out_rrp_data = bank[0].port_readdata[4];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[4];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port4bank1
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [3:0] icm_in_arb_address [1];
            logic [63:0] icm_in_arb_writedata [1];
            logic [7:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [63:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [3:0] icm_out_arb_address;
            logic [63:0] icm_out_arb_writedata;
            logic [7:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [63:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[5].b_arb_request[1];
            assign icm_in_arb_enable[0] = router[5].b_arb_enable[1];
            assign icm_in_arb_read[0] = router[5].b_arb_read[1];
            assign icm_in_arb_write[0] = router[5].b_arb_write[1];
            assign icm_in_arb_burstcount[0] = router[5].b_arb_burstcount[1];
            assign icm_in_arb_address[0] = router[5].b_arb_address[1];
            assign icm_in_arb_writedata[0] = router[5].b_arb_writedata[1];
            assign icm_in_arb_byteenable[0] = router[5].b_arb_byteenable[1];
            assign router[5].b_arb_stall[1] = icm_in_arb_stall[0];
            assign router[5].b_wrp_ack[1] = icm_in_wrp_ack[0];
            assign router[5].b_rrp_datavalid[1] = icm_in_rrp_datavalid[0];
            assign router[5].b_rrp_data[1] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_1
            Krnl_GA_system_interconnect_1 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[1].port_enable[4] = icm_out_arb_enable;
            assign bank[1].port_read[4] = icm_out_arb_read;
            assign bank[1].port_write[4] = icm_out_arb_write;
            assign bank[1].port_address[4] = icm_out_arb_address;
            assign bank[1].port_writedata[4] = icm_out_arb_writedata;
            assign bank[1].port_byteenable[4] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[1].port_waitrequest[4];
            assign icm_out_rrp_data = bank[1].port_readdata[4];
            assign icm_out_rrp_datavalid = bank[1].port_readdatavalid[4];
            assign icm_out_wrp_ack = 'b0;
         end

      end

      assign lmem_invalid_aspaces[0] = |invalid_access_grps;
   end
   endgenerate

   generate
   begin:local_mem_system_aspace7
      logic local_icm_arb_request [1][3];
      logic local_icm_arb_enable [1][3];
      logic local_icm_arb_read [1][3];
      logic local_icm_arb_write [1][3];
      logic local_icm_arb_burstcount [1][3];
      logic [6:0] local_icm_arb_address [1][3];
      logic [31:0] local_icm_arb_writedata [1][3];
      logic [3:0] local_icm_arb_byteenable [1][3];
      logic local_icm_arb_stall [1][3];
      logic local_icm_wrp_ack [1][3];
      logic local_icm_rrp_datavalid [1][3];
      logic [31:0] local_icm_rrp_data [1][3];

      for( __j = 0; __j < 1; __j = __j + 1 )
      begin:local_mem_group
         for( __k = 0; __k < 3; __k = __k + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace7_enable[__j][__k]),
               .avm_read(local_avm_aspace7_read[__j][__k]),
               .avm_write(local_avm_aspace7_write[__j][__k]),
               .avm_burstcount(local_avm_aspace7_burstcount[__j][__k]),
               .avm_address(local_avm_aspace7_address[__j][__k]),
               .avm_writedata(local_avm_aspace7_writedata[__j][__k]),
               .avm_byteenable(local_avm_aspace7_byteenable[__j][__k]),
               .avm_waitrequest(local_avm_aspace7_waitrequest[__j][__k]),
               .avm_readdata(local_avm_aspace7_readdata[__j][__k]),
               .avm_readdatavalid(local_avm_aspace7_readdatavalid[__j][__k]),
               .avm_writeack(local_avm_aspace7_writeack[__j][__k]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__j][__k]),
               .ic_arb_enable(local_icm_arb_enable[__j][__k]),
               .ic_arb_read(local_icm_arb_read[__j][__k]),
               .ic_arb_write(local_icm_arb_write[__j][__k]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__j][__k]),
               .ic_arb_address(local_icm_arb_address[__j][__k]),
               .ic_arb_writedata(local_icm_arb_writedata[__j][__k]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__j][__k]),
               .ic_arb_stall(local_icm_arb_stall[__j][__k]),
               .ic_wrp_ack(local_icm_wrp_ack[__j][__k]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__j][__k]),
               .ic_rrp_data(local_icm_rrp_data[__j][__k])
            );

         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:bank
            logic port_enable [1:4];
            logic port_read [1:4];
            logic port_write [1:4];
            logic [6:0] port_address [1:4];
            logic [31:0] port_writedata [1:4];
            logic [3:0] port_byteenable [1:4];
            logic port_waitrequest [1:4];
            logic [31:0] port_readdata [1:4];
            logic port_readdatavalid [1:4];

            // INST mem0 of acl_mem2x
            acl_mem2x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("BIDIR_DUAL_PORT"),
               .PREFERRED_WIDTH(160),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .clk2x(clock2x),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2]),
               // AVS avs_port3
               .avs_port3_enable(port_enable[3]),
               .avs_port3_read(port_read[3]),
               .avs_port3_write(port_write[3]),
               .avs_port3_address(port_address[3]),
               .avs_port3_writedata(port_writedata[3]),
               .avs_port3_byteenable(port_byteenable[3]),
               .avs_port3_waitrequest(port_waitrequest[3]),
               .avs_port3_readdata(port_readdata[3]),
               .avs_port3_readdatavalid(port_readdatavalid[3]),
               // AVS avs_port4
               .avs_port4_enable(port_enable[4]),
               .avs_port4_read(port_read[4]),
               .avs_port4_write(port_write[4]),
               .avs_port4_address(port_address[4]),
               .avs_port4_writedata(port_writedata[4]),
               .avs_port4_byteenable(port_byteenable[4]),
               .avs_port4_waitrequest(port_waitrequest[4]),
               .avs_port4_readdata(port_readdata[4]),
               .avs_port4_readdatavalid(port_readdatavalid[4])
            );

         end

         for( __k = 0; __k < 3; __k = __k + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__j][__k]),
               .m_arb_enable(local_icm_arb_enable[__j][__k]),
               .m_arb_read(local_icm_arb_read[__j][__k]),
               .m_arb_write(local_icm_arb_write[__j][__k]),
               .m_arb_burstcount(local_icm_arb_burstcount[__j][__k]),
               .m_arb_address(local_icm_arb_address[__j][__k]),
               .m_arb_writedata(local_icm_arb_writedata[__j][__k]),
               .m_arb_byteenable(local_icm_arb_byteenable[__j][__k]),
               .m_arb_stall(local_icm_arb_stall[__j][__k]),
               .m_wrp_ack(local_icm_wrp_ack[__j][__k]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__j][__k]),
               .m_rrp_data(local_icm_rrp_data[__j][__k]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_3
            Krnl_GA_system_interconnect_3 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port3bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[2].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[2].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[2].b_arb_read[0];
            assign icm_in_arb_write[0] = router[2].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[2].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[2].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[2].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[2].b_arb_byteenable[0];
            assign router[2].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[2].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[2].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[2].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[3] = icm_out_arb_enable;
            assign bank[0].port_read[3] = icm_out_arb_read;
            assign bank[0].port_write[3] = icm_out_arb_write;
            assign bank[0].port_address[3] = icm_out_arb_address;
            assign bank[0].port_writedata[3] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[3];
            assign icm_out_rrp_data = bank[0].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port4bank0
            assign bank[0].port_enable[4] = '0;
            assign bank[0].port_read[4] = '0;
            assign bank[0].port_write[4] = '0;
            assign bank[0].port_address[4] = '0;
            assign bank[0].port_writedata[4] = '0;
            assign bank[0].port_byteenable[4] = '0;
         end

      end

   end
   endgenerate

   generate
   begin:local_mem_system_aspace8
      logic local_icm_arb_request [1][3];
      logic local_icm_arb_enable [1][3];
      logic local_icm_arb_read [1][3];
      logic local_icm_arb_write [1][3];
      logic local_icm_arb_burstcount [1][3];
      logic [6:0] local_icm_arb_address [1][3];
      logic [31:0] local_icm_arb_writedata [1][3];
      logic [3:0] local_icm_arb_byteenable [1][3];
      logic local_icm_arb_stall [1][3];
      logic local_icm_wrp_ack [1][3];
      logic local_icm_rrp_datavalid [1][3];
      logic [31:0] local_icm_rrp_data [1][3];

      for( __k = 0; __k < 1; __k = __k + 1 )
      begin:local_mem_group
         for( __l = 0; __l < 3; __l = __l + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace8_enable[__k][__l]),
               .avm_read(local_avm_aspace8_read[__k][__l]),
               .avm_write(local_avm_aspace8_write[__k][__l]),
               .avm_burstcount(local_avm_aspace8_burstcount[__k][__l]),
               .avm_address(local_avm_aspace8_address[__k][__l]),
               .avm_writedata(local_avm_aspace8_writedata[__k][__l]),
               .avm_byteenable(local_avm_aspace8_byteenable[__k][__l]),
               .avm_waitrequest(local_avm_aspace8_waitrequest[__k][__l]),
               .avm_readdata(local_avm_aspace8_readdata[__k][__l]),
               .avm_readdatavalid(local_avm_aspace8_readdatavalid[__k][__l]),
               .avm_writeack(local_avm_aspace8_writeack[__k][__l]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__k][__l]),
               .ic_arb_enable(local_icm_arb_enable[__k][__l]),
               .ic_arb_read(local_icm_arb_read[__k][__l]),
               .ic_arb_write(local_icm_arb_write[__k][__l]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__k][__l]),
               .ic_arb_address(local_icm_arb_address[__k][__l]),
               .ic_arb_writedata(local_icm_arb_writedata[__k][__l]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__k][__l]),
               .ic_arb_stall(local_icm_arb_stall[__k][__l]),
               .ic_wrp_ack(local_icm_wrp_ack[__k][__l]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__k][__l]),
               .ic_rrp_data(local_icm_rrp_data[__k][__l])
            );

         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:bank
            logic port_enable [1:4];
            logic port_read [1:4];
            logic port_write [1:4];
            logic [6:0] port_address [1:4];
            logic [31:0] port_writedata [1:4];
            logic [3:0] port_byteenable [1:4];
            logic port_waitrequest [1:4];
            logic [31:0] port_readdata [1:4];
            logic port_readdatavalid [1:4];

            // INST mem0 of acl_mem2x
            acl_mem2x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("BIDIR_DUAL_PORT"),
               .PREFERRED_WIDTH(160),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .clk2x(clock2x),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2]),
               // AVS avs_port3
               .avs_port3_enable(port_enable[3]),
               .avs_port3_read(port_read[3]),
               .avs_port3_write(port_write[3]),
               .avs_port3_address(port_address[3]),
               .avs_port3_writedata(port_writedata[3]),
               .avs_port3_byteenable(port_byteenable[3]),
               .avs_port3_waitrequest(port_waitrequest[3]),
               .avs_port3_readdata(port_readdata[3]),
               .avs_port3_readdatavalid(port_readdatavalid[3]),
               // AVS avs_port4
               .avs_port4_enable(port_enable[4]),
               .avs_port4_read(port_read[4]),
               .avs_port4_write(port_write[4]),
               .avs_port4_address(port_address[4]),
               .avs_port4_writedata(port_writedata[4]),
               .avs_port4_byteenable(port_byteenable[4]),
               .avs_port4_waitrequest(port_waitrequest[4]),
               .avs_port4_readdata(port_readdata[4]),
               .avs_port4_readdatavalid(port_readdatavalid[4])
            );

         end

         for( __l = 0; __l < 3; __l = __l + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__k][__l]),
               .m_arb_enable(local_icm_arb_enable[__k][__l]),
               .m_arb_read(local_icm_arb_read[__k][__l]),
               .m_arb_write(local_icm_arb_write[__k][__l]),
               .m_arb_burstcount(local_icm_arb_burstcount[__k][__l]),
               .m_arb_address(local_icm_arb_address[__k][__l]),
               .m_arb_writedata(local_icm_arb_writedata[__k][__l]),
               .m_arb_byteenable(local_icm_arb_byteenable[__k][__l]),
               .m_arb_stall(local_icm_arb_stall[__k][__l]),
               .m_wrp_ack(local_icm_wrp_ack[__k][__l]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__k][__l]),
               .m_rrp_data(local_icm_rrp_data[__k][__l]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_3
            Krnl_GA_system_interconnect_3 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port3bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[2].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[2].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[2].b_arb_read[0];
            assign icm_in_arb_write[0] = router[2].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[2].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[2].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[2].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[2].b_arb_byteenable[0];
            assign router[2].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[2].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[2].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[2].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[3] = icm_out_arb_enable;
            assign bank[0].port_read[3] = icm_out_arb_read;
            assign bank[0].port_write[3] = icm_out_arb_write;
            assign bank[0].port_address[3] = icm_out_arb_address;
            assign bank[0].port_writedata[3] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[3];
            assign icm_out_rrp_data = bank[0].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port4bank0
            assign bank[0].port_enable[4] = '0;
            assign bank[0].port_read[4] = '0;
            assign bank[0].port_write[4] = '0;
            assign bank[0].port_address[4] = '0;
            assign bank[0].port_writedata[4] = '0;
            assign bank[0].port_byteenable[4] = '0;
         end

      end

   end
   endgenerate

   generate
   begin:local_mem_system_aspace9
      logic local_icm_arb_request [1][3];
      logic local_icm_arb_enable [1][3];
      logic local_icm_arb_read [1][3];
      logic local_icm_arb_write [1][3];
      logic local_icm_arb_burstcount [1][3];
      logic [6:0] local_icm_arb_address [1][3];
      logic [31:0] local_icm_arb_writedata [1][3];
      logic [3:0] local_icm_arb_byteenable [1][3];
      logic local_icm_arb_stall [1][3];
      logic local_icm_wrp_ack [1][3];
      logic local_icm_rrp_datavalid [1][3];
      logic [31:0] local_icm_rrp_data [1][3];

      for( __l = 0; __l < 1; __l = __l + 1 )
      begin:local_mem_group
         for( __m = 0; __m < 3; __m = __m + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace9_enable[__l][__m]),
               .avm_read(local_avm_aspace9_read[__l][__m]),
               .avm_write(local_avm_aspace9_write[__l][__m]),
               .avm_burstcount(local_avm_aspace9_burstcount[__l][__m]),
               .avm_address(local_avm_aspace9_address[__l][__m]),
               .avm_writedata(local_avm_aspace9_writedata[__l][__m]),
               .avm_byteenable(local_avm_aspace9_byteenable[__l][__m]),
               .avm_waitrequest(local_avm_aspace9_waitrequest[__l][__m]),
               .avm_readdata(local_avm_aspace9_readdata[__l][__m]),
               .avm_readdatavalid(local_avm_aspace9_readdatavalid[__l][__m]),
               .avm_writeack(local_avm_aspace9_writeack[__l][__m]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__l][__m]),
               .ic_arb_enable(local_icm_arb_enable[__l][__m]),
               .ic_arb_read(local_icm_arb_read[__l][__m]),
               .ic_arb_write(local_icm_arb_write[__l][__m]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__l][__m]),
               .ic_arb_address(local_icm_arb_address[__l][__m]),
               .ic_arb_writedata(local_icm_arb_writedata[__l][__m]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__l][__m]),
               .ic_arb_stall(local_icm_arb_stall[__l][__m]),
               .ic_wrp_ack(local_icm_wrp_ack[__l][__m]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__l][__m]),
               .ic_rrp_data(local_icm_rrp_data[__l][__m])
            );

         end

         for( __m = 0; __m < 1; __m = __m + 1 )
         begin:bank
            logic port_enable [1:4];
            logic port_read [1:4];
            logic port_write [1:4];
            logic [6:0] port_address [1:4];
            logic [31:0] port_writedata [1:4];
            logic [3:0] port_byteenable [1:4];
            logic port_waitrequest [1:4];
            logic [31:0] port_readdata [1:4];
            logic port_readdatavalid [1:4];

            // INST mem0 of acl_mem2x
            acl_mem2x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("BIDIR_DUAL_PORT"),
               .PREFERRED_WIDTH(160),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .clk2x(clock2x),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2]),
               // AVS avs_port3
               .avs_port3_enable(port_enable[3]),
               .avs_port3_read(port_read[3]),
               .avs_port3_write(port_write[3]),
               .avs_port3_address(port_address[3]),
               .avs_port3_writedata(port_writedata[3]),
               .avs_port3_byteenable(port_byteenable[3]),
               .avs_port3_waitrequest(port_waitrequest[3]),
               .avs_port3_readdata(port_readdata[3]),
               .avs_port3_readdatavalid(port_readdatavalid[3]),
               // AVS avs_port4
               .avs_port4_enable(port_enable[4]),
               .avs_port4_read(port_read[4]),
               .avs_port4_write(port_write[4]),
               .avs_port4_address(port_address[4]),
               .avs_port4_writedata(port_writedata[4]),
               .avs_port4_byteenable(port_byteenable[4]),
               .avs_port4_waitrequest(port_waitrequest[4]),
               .avs_port4_readdata(port_readdata[4]),
               .avs_port4_readdatavalid(port_readdatavalid[4])
            );

         end

         for( __m = 0; __m < 3; __m = __m + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__l][__m]),
               .m_arb_enable(local_icm_arb_enable[__l][__m]),
               .m_arb_read(local_icm_arb_read[__l][__m]),
               .m_arb_write(local_icm_arb_write[__l][__m]),
               .m_arb_burstcount(local_icm_arb_burstcount[__l][__m]),
               .m_arb_address(local_icm_arb_address[__l][__m]),
               .m_arb_writedata(local_icm_arb_writedata[__l][__m]),
               .m_arb_byteenable(local_icm_arb_byteenable[__l][__m]),
               .m_arb_stall(local_icm_arb_stall[__l][__m]),
               .m_wrp_ack(local_icm_wrp_ack[__l][__m]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__l][__m]),
               .m_rrp_data(local_icm_rrp_data[__l][__m]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __m = 0; __m < 1; __m = __m + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __m = 0; __m < 1; __m = __m + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_3
            Krnl_GA_system_interconnect_3 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __m = 0; __m < 1; __m = __m + 1 )
         begin:port3bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[2].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[2].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[2].b_arb_read[0];
            assign icm_in_arb_write[0] = router[2].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[2].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[2].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[2].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[2].b_arb_byteenable[0];
            assign router[2].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[2].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[2].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[2].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[3] = icm_out_arb_enable;
            assign bank[0].port_read[3] = icm_out_arb_read;
            assign bank[0].port_write[3] = icm_out_arb_write;
            assign bank[0].port_address[3] = icm_out_arb_address;
            assign bank[0].port_writedata[3] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[3];
            assign icm_out_rrp_data = bank[0].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __m = 0; __m < 1; __m = __m + 1 )
         begin:port4bank0
            assign bank[0].port_enable[4] = '0;
            assign bank[0].port_read[4] = '0;
            assign bank[0].port_write[4] = '0;
            assign bank[0].port_address[4] = '0;
            assign bank[0].port_writedata[4] = '0;
            assign bank[0].port_byteenable[4] = '0;
         end

      end

   end
   endgenerate

endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_top_wrapper_0
/////////////////////////////////////////////////////////////////
module Krnl_GA_top_wrapper_0
(
   input logic start,
   input logic [671:0] kernel_arguments,
   input logic [31:0] work_dim,
   input logic [31:0] global_offset [2:0],
   output logic kernel_valid_out,
   output logic has_a_write_pending,
   output logic has_a_lsu_active,
   input logic [31:0] global_id [2:0],
   input logic [31:0] local_id [2:0],
   input logic [31:0] group_id [2:0],
   input logic [31:0] global_size [2:0],
   input logic [31:0] local_size [2:0],
   input logic [31:0] num_groups [2:0],
   input logic [31:0] workgroup_size,
   output logic kernel_stall_out,
   input logic kernel_valid_in,
   input logic clock,
   input logic resetn,
   input logic clock2x,
   // AVM avm_local_bb0_st__inst0
   output logic avm_local_bb0_st__inst0_enable,
   output logic avm_local_bb0_st__inst0_read,
   output logic avm_local_bb0_st__inst0_write,
   output logic [4:0] avm_local_bb0_st__inst0_burstcount,
   output logic [30:0] avm_local_bb0_st__inst0_address,
   output logic [511:0] avm_local_bb0_st__inst0_writedata,
   output logic [63:0] avm_local_bb0_st__inst0_byteenable,
   input logic avm_local_bb0_st__inst0_waitrequest,
   input logic [511:0] avm_local_bb0_st__inst0_readdata,
   input logic avm_local_bb0_st__inst0_readdatavalid,
   input logic avm_local_bb0_st__inst0_writeack,
   // AVM avm_local_bb2_st__inst0
   output logic avm_local_bb2_st__inst0_enable,
   output logic avm_local_bb2_st__inst0_read,
   output logic avm_local_bb2_st__inst0_write,
   output logic [4:0] avm_local_bb2_st__inst0_burstcount,
   output logic [30:0] avm_local_bb2_st__inst0_address,
   output logic [511:0] avm_local_bb2_st__inst0_writedata,
   output logic [63:0] avm_local_bb2_st__inst0_byteenable,
   input logic avm_local_bb2_st__inst0_waitrequest,
   input logic [511:0] avm_local_bb2_st__inst0_readdata,
   input logic avm_local_bb2_st__inst0_readdatavalid,
   input logic avm_local_bb2_st__inst0_writeack,
   // AVM avm_local_bb2_st__u1_inst0
   output logic avm_local_bb2_st__u1_inst0_enable,
   output logic avm_local_bb2_st__u1_inst0_read,
   output logic avm_local_bb2_st__u1_inst0_write,
   output logic [4:0] avm_local_bb2_st__u1_inst0_burstcount,
   output logic [30:0] avm_local_bb2_st__u1_inst0_address,
   output logic [511:0] avm_local_bb2_st__u1_inst0_writedata,
   output logic [63:0] avm_local_bb2_st__u1_inst0_byteenable,
   input logic avm_local_bb2_st__u1_inst0_waitrequest,
   input logic [511:0] avm_local_bb2_st__u1_inst0_readdata,
   input logic avm_local_bb2_st__u1_inst0_readdatavalid,
   input logic avm_local_bb2_st__u1_inst0_writeack,
   // AVST avm_channel_id_chan_GA2Conf_genotype_write
   output logic avm_channel_id_chan_GA2Conf_genotype_write_valid,
   input logic avm_channel_id_chan_GA2Conf_genotype_write_ready,
   output logic [31:0] avm_channel_id_chan_GA2Conf_genotype_write_data,
   input logic avm_channel_id_chan_GA2Conf_genotype_write_almostfull,
   // AVST avm_channel_id_chan_GA2Conf_active_write
   output logic avm_channel_id_chan_GA2Conf_active_write_valid,
   input logic avm_channel_id_chan_GA2Conf_active_write_ready,
   output logic [7:0] avm_channel_id_chan_GA2Conf_active_write_data,
   input logic avm_channel_id_chan_GA2Conf_active_write_almostfull,
   // AVST avm_channel_id_chan_GA2Conf_cnt_write
   output logic avm_channel_id_chan_GA2Conf_cnt_write_valid,
   input logic avm_channel_id_chan_GA2Conf_cnt_write_ready,
   output logic [31:0] avm_channel_id_chan_GA2Conf_cnt_write_data,
   input logic avm_channel_id_chan_GA2Conf_cnt_write_almostfull,
   // AVST avm_channel_id_chan_GA2Conf_mode_write
   output logic avm_channel_id_chan_GA2Conf_mode_write_valid,
   input logic avm_channel_id_chan_GA2Conf_mode_write_ready,
   output logic [7:0] avm_channel_id_chan_GA2Conf_mode_write_data,
   input logic avm_channel_id_chan_GA2Conf_mode_write_almostfull,
   // AVM p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid
);
   logic lmem_invalid_single_bit;

   // INST kernel of Krnl_GA_function_wrapper
   Krnl_GA_function_wrapper kernel
   (
      .local_router_hang(lmem_invalid_single_bit),
      .start(start),
      .kernel_arguments(kernel_arguments),
      .work_dim(work_dim),
      .global_offset_0(global_offset[0]),
      .global_offset_1(global_offset[1]),
      .global_offset_2(global_offset[2]),
      .kernel_valid_out(kernel_valid_out),
      .has_a_write_pending(has_a_write_pending),
      .has_a_lsu_active(has_a_lsu_active),
      .global_id_0(global_id[0]),
      .global_id_1(global_id[1]),
      .global_id_2(global_id[2]),
      .local_id_0(local_id[0]),
      .local_id_1(local_id[1]),
      .local_id_2(local_id[2]),
      .group_id_0(group_id[0]),
      .group_id_1(group_id[1]),
      .group_id_2(group_id[2]),
      .global_size_0(global_size[0]),
      .global_size_1(global_size[1]),
      .global_size_2(global_size[2]),
      .local_size_0(local_size[0]),
      .local_size_1(local_size[1]),
      .local_size_2(local_size[2]),
      .num_groups_0(num_groups[0]),
      .num_groups_1(num_groups[1]),
      .num_groups_2(num_groups[2]),
      .workgroup_size(workgroup_size),
      .kernel_stall_out(kernel_stall_out),
      .kernel_valid_in(kernel_valid_in),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_st__inst0
      .avm_local_bb0_st__inst0_enable(avm_local_bb0_st__inst0_enable),
      .avm_local_bb0_st__inst0_read(avm_local_bb0_st__inst0_read),
      .avm_local_bb0_st__inst0_write(avm_local_bb0_st__inst0_write),
      .avm_local_bb0_st__inst0_burstcount(avm_local_bb0_st__inst0_burstcount),
      .avm_local_bb0_st__inst0_address(avm_local_bb0_st__inst0_address),
      .avm_local_bb0_st__inst0_writedata(avm_local_bb0_st__inst0_writedata),
      .avm_local_bb0_st__inst0_byteenable(avm_local_bb0_st__inst0_byteenable),
      .avm_local_bb0_st__inst0_waitrequest(avm_local_bb0_st__inst0_waitrequest),
      .avm_local_bb0_st__inst0_readdata(avm_local_bb0_st__inst0_readdata),
      .avm_local_bb0_st__inst0_readdatavalid(avm_local_bb0_st__inst0_readdatavalid),
      .avm_local_bb0_st__inst0_writeack(avm_local_bb0_st__inst0_writeack),
      // AVM avm_local_bb2_st__inst0
      .avm_local_bb2_st__inst0_enable(avm_local_bb2_st__inst0_enable),
      .avm_local_bb2_st__inst0_read(avm_local_bb2_st__inst0_read),
      .avm_local_bb2_st__inst0_write(avm_local_bb2_st__inst0_write),
      .avm_local_bb2_st__inst0_burstcount(avm_local_bb2_st__inst0_burstcount),
      .avm_local_bb2_st__inst0_address(avm_local_bb2_st__inst0_address),
      .avm_local_bb2_st__inst0_writedata(avm_local_bb2_st__inst0_writedata),
      .avm_local_bb2_st__inst0_byteenable(avm_local_bb2_st__inst0_byteenable),
      .avm_local_bb2_st__inst0_waitrequest(avm_local_bb2_st__inst0_waitrequest),
      .avm_local_bb2_st__inst0_readdata(avm_local_bb2_st__inst0_readdata),
      .avm_local_bb2_st__inst0_readdatavalid(avm_local_bb2_st__inst0_readdatavalid),
      .avm_local_bb2_st__inst0_writeack(avm_local_bb2_st__inst0_writeack),
      // AVM avm_local_bb2_st__u1_inst0
      .avm_local_bb2_st__u1_inst0_enable(avm_local_bb2_st__u1_inst0_enable),
      .avm_local_bb2_st__u1_inst0_read(avm_local_bb2_st__u1_inst0_read),
      .avm_local_bb2_st__u1_inst0_write(avm_local_bb2_st__u1_inst0_write),
      .avm_local_bb2_st__u1_inst0_burstcount(avm_local_bb2_st__u1_inst0_burstcount),
      .avm_local_bb2_st__u1_inst0_address(avm_local_bb2_st__u1_inst0_address),
      .avm_local_bb2_st__u1_inst0_writedata(avm_local_bb2_st__u1_inst0_writedata),
      .avm_local_bb2_st__u1_inst0_byteenable(avm_local_bb2_st__u1_inst0_byteenable),
      .avm_local_bb2_st__u1_inst0_waitrequest(avm_local_bb2_st__u1_inst0_waitrequest),
      .avm_local_bb2_st__u1_inst0_readdata(avm_local_bb2_st__u1_inst0_readdata),
      .avm_local_bb2_st__u1_inst0_readdatavalid(avm_local_bb2_st__u1_inst0_readdatavalid),
      .avm_local_bb2_st__u1_inst0_writeack(avm_local_bb2_st__u1_inst0_writeack),
      // AVST avst_local_bb1__chan_GA2Conf_genotype_inst0
      .avst_local_bb1__chan_GA2Conf_genotype_inst0_valid(avm_channel_id_chan_GA2Conf_genotype_write_valid),
      .avst_local_bb1__chan_GA2Conf_genotype_inst0_ready(avm_channel_id_chan_GA2Conf_genotype_write_ready),
      .avst_local_bb1__chan_GA2Conf_genotype_inst0_data(avm_channel_id_chan_GA2Conf_genotype_write_data),
      .avst_local_bb1__chan_GA2Conf_genotype_inst0_almostfull(avm_channel_id_chan_GA2Conf_genotype_write_almostfull),
      // AVST avst_local_bb2__chan_GA2Conf_active_inst0
      .avst_local_bb2__chan_GA2Conf_active_inst0_valid(avm_channel_id_chan_GA2Conf_active_write_valid),
      .avst_local_bb2__chan_GA2Conf_active_inst0_ready(avm_channel_id_chan_GA2Conf_active_write_ready),
      .avst_local_bb2__chan_GA2Conf_active_inst0_data(avm_channel_id_chan_GA2Conf_active_write_data),
      .avst_local_bb2__chan_GA2Conf_active_inst0_almostfull(avm_channel_id_chan_GA2Conf_active_write_almostfull),
      // AVST avst_local_bb2__chan_GA2Conf_cnt_inst0
      .avst_local_bb2__chan_GA2Conf_cnt_inst0_valid(avm_channel_id_chan_GA2Conf_cnt_write_valid),
      .avst_local_bb2__chan_GA2Conf_cnt_inst0_ready(avm_channel_id_chan_GA2Conf_cnt_write_ready),
      .avst_local_bb2__chan_GA2Conf_cnt_inst0_data(avm_channel_id_chan_GA2Conf_cnt_write_data),
      .avst_local_bb2__chan_GA2Conf_cnt_inst0_almostfull(avm_channel_id_chan_GA2Conf_cnt_write_almostfull),
      // AVST avst_local_bb2__chan_GA2Conf_mode_inst0
      .avst_local_bb2__chan_GA2Conf_mode_inst0_valid(avm_channel_id_chan_GA2Conf_mode_write_valid),
      .avst_local_bb2__chan_GA2Conf_mode_inst0_ready(avm_channel_id_chan_GA2Conf_mode_write_ready),
      .avst_local_bb2__chan_GA2Conf_mode_inst0_data(avm_channel_id_chan_GA2Conf_mode_write_data),
      .avst_local_bb2__chan_GA2Conf_mode_inst0_almostfull(avm_channel_id_chan_GA2Conf_mode_write_almostfull),
      // AVM p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_read(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_write(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_address(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb0_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid)
   );

   assign lmem_invalid_single_bit = 'b0;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_InterE_top_wrapper_0
/////////////////////////////////////////////////////////////////
module Krnl_InterE_top_wrapper_0
(
   input logic start,
   input logic [287:0] kernel_arguments,
   input logic [31:0] work_dim,
   input logic [31:0] global_offset [2:0],
   output logic kernel_valid_out,
   output logic has_a_write_pending,
   output logic has_a_lsu_active,
   input logic [31:0] global_id [2:0],
   input logic [31:0] local_id [2:0],
   input logic [31:0] group_id [2:0],
   input logic [31:0] global_size [2:0],
   input logic [31:0] local_size [2:0],
   input logic [31:0] num_groups [2:0],
   input logic [31:0] workgroup_size,
   output logic kernel_stall_out,
   input logic kernel_valid_in,
   input logic clock,
   input logic resetn,
   input logic clock2x,
   // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable,
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read,
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write,
   output logic [4:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount,
   output logic [30:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address,
   output logic [511:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata,
   output logic [63:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest,
   input logic [511:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack,
   // AVM avm_local_bb3_st__inst0
   output logic avm_local_bb3_st__inst0_enable,
   output logic avm_local_bb3_st__inst0_read,
   output logic avm_local_bb3_st__inst0_write,
   output logic [4:0] avm_local_bb3_st__inst0_burstcount,
   output logic [30:0] avm_local_bb3_st__inst0_address,
   output logic [511:0] avm_local_bb3_st__inst0_writedata,
   output logic [63:0] avm_local_bb3_st__inst0_byteenable,
   input logic avm_local_bb3_st__inst0_waitrequest,
   input logic [511:0] avm_local_bb3_st__inst0_readdata,
   input logic avm_local_bb3_st__inst0_readdatavalid,
   input logic avm_local_bb3_st__inst0_writeack,
   // AVM avm_local_bb4_ld__inst0
   output logic avm_local_bb4_ld__inst0_enable,
   output logic avm_local_bb4_ld__inst0_read,
   output logic avm_local_bb4_ld__inst0_write,
   output logic [4:0] avm_local_bb4_ld__inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__inst0_address,
   output logic [511:0] avm_local_bb4_ld__inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__inst0_byteenable,
   input logic avm_local_bb4_ld__inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__inst0_readdata,
   input logic avm_local_bb4_ld__inst0_readdatavalid,
   input logic avm_local_bb4_ld__inst0_writeack,
   // AVM avm_local_bb4_ld__u11_inst0
   output logic avm_local_bb4_ld__u11_inst0_enable,
   output logic avm_local_bb4_ld__u11_inst0_read,
   output logic avm_local_bb4_ld__u11_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u11_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u11_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u11_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u11_inst0_byteenable,
   input logic avm_local_bb4_ld__u11_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u11_inst0_readdata,
   input logic avm_local_bb4_ld__u11_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u11_inst0_writeack,
   // AVM avm_local_bb4_ld__u12_inst0
   output logic avm_local_bb4_ld__u12_inst0_enable,
   output logic avm_local_bb4_ld__u12_inst0_read,
   output logic avm_local_bb4_ld__u12_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u12_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u12_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u12_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u12_inst0_byteenable,
   input logic avm_local_bb4_ld__u12_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u12_inst0_readdata,
   input logic avm_local_bb4_ld__u12_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u12_inst0_writeack,
   // AVM avm_local_bb4_ld__u13_inst0
   output logic avm_local_bb4_ld__u13_inst0_enable,
   output logic avm_local_bb4_ld__u13_inst0_read,
   output logic avm_local_bb4_ld__u13_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u13_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u13_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u13_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u13_inst0_byteenable,
   input logic avm_local_bb4_ld__u13_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u13_inst0_readdata,
   input logic avm_local_bb4_ld__u13_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u13_inst0_writeack,
   // AVM avm_local_bb4_ld__u14_inst0
   output logic avm_local_bb4_ld__u14_inst0_enable,
   output logic avm_local_bb4_ld__u14_inst0_read,
   output logic avm_local_bb4_ld__u14_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u14_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u14_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u14_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u14_inst0_byteenable,
   input logic avm_local_bb4_ld__u14_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u14_inst0_readdata,
   input logic avm_local_bb4_ld__u14_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u14_inst0_writeack,
   // AVM avm_local_bb4_ld__u15_inst0
   output logic avm_local_bb4_ld__u15_inst0_enable,
   output logic avm_local_bb4_ld__u15_inst0_read,
   output logic avm_local_bb4_ld__u15_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u15_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u15_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u15_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u15_inst0_byteenable,
   input logic avm_local_bb4_ld__u15_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u15_inst0_readdata,
   input logic avm_local_bb4_ld__u15_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u15_inst0_writeack,
   // AVM avm_local_bb4_ld__u16_inst0
   output logic avm_local_bb4_ld__u16_inst0_enable,
   output logic avm_local_bb4_ld__u16_inst0_read,
   output logic avm_local_bb4_ld__u16_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u16_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u16_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u16_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u16_inst0_byteenable,
   input logic avm_local_bb4_ld__u16_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u16_inst0_readdata,
   input logic avm_local_bb4_ld__u16_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u16_inst0_writeack,
   // AVM avm_local_bb4_ld__u17_inst0
   output logic avm_local_bb4_ld__u17_inst0_enable,
   output logic avm_local_bb4_ld__u17_inst0_read,
   output logic avm_local_bb4_ld__u17_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u17_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u17_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u17_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u17_inst0_byteenable,
   input logic avm_local_bb4_ld__u17_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u17_inst0_readdata,
   input logic avm_local_bb4_ld__u17_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u17_inst0_writeack,
   // AVM avm_local_bb4_ld__u18_inst0
   output logic avm_local_bb4_ld__u18_inst0_enable,
   output logic avm_local_bb4_ld__u18_inst0_read,
   output logic avm_local_bb4_ld__u18_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u18_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u18_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u18_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u18_inst0_byteenable,
   input logic avm_local_bb4_ld__u18_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u18_inst0_readdata,
   input logic avm_local_bb4_ld__u18_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u18_inst0_writeack,
   // AVM avm_local_bb4_ld__u19_inst0
   output logic avm_local_bb4_ld__u19_inst0_enable,
   output logic avm_local_bb4_ld__u19_inst0_read,
   output logic avm_local_bb4_ld__u19_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u19_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u19_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u19_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u19_inst0_byteenable,
   input logic avm_local_bb4_ld__u19_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u19_inst0_readdata,
   input logic avm_local_bb4_ld__u19_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u19_inst0_writeack,
   // AVM avm_local_bb4_ld__u20_inst0
   output logic avm_local_bb4_ld__u20_inst0_enable,
   output logic avm_local_bb4_ld__u20_inst0_read,
   output logic avm_local_bb4_ld__u20_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u20_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u20_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u20_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u20_inst0_byteenable,
   input logic avm_local_bb4_ld__u20_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u20_inst0_readdata,
   input logic avm_local_bb4_ld__u20_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u20_inst0_writeack,
   // AVM avm_local_bb4_ld__u21_inst0
   output logic avm_local_bb4_ld__u21_inst0_enable,
   output logic avm_local_bb4_ld__u21_inst0_read,
   output logic avm_local_bb4_ld__u21_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u21_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u21_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u21_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u21_inst0_byteenable,
   input logic avm_local_bb4_ld__u21_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u21_inst0_readdata,
   input logic avm_local_bb4_ld__u21_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u21_inst0_writeack,
   // AVM avm_local_bb4_ld__u22_inst0
   output logic avm_local_bb4_ld__u22_inst0_enable,
   output logic avm_local_bb4_ld__u22_inst0_read,
   output logic avm_local_bb4_ld__u22_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u22_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u22_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u22_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u22_inst0_byteenable,
   input logic avm_local_bb4_ld__u22_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u22_inst0_readdata,
   input logic avm_local_bb4_ld__u22_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u22_inst0_writeack,
   // AVM avm_local_bb4_ld__u23_inst0
   output logic avm_local_bb4_ld__u23_inst0_enable,
   output logic avm_local_bb4_ld__u23_inst0_read,
   output logic avm_local_bb4_ld__u23_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u23_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u23_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u23_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u23_inst0_byteenable,
   input logic avm_local_bb4_ld__u23_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u23_inst0_readdata,
   input logic avm_local_bb4_ld__u23_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u23_inst0_writeack,
   // AVM avm_local_bb4_ld__u24_inst0
   output logic avm_local_bb4_ld__u24_inst0_enable,
   output logic avm_local_bb4_ld__u24_inst0_read,
   output logic avm_local_bb4_ld__u24_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u24_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u24_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u24_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u24_inst0_byteenable,
   input logic avm_local_bb4_ld__u24_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u24_inst0_readdata,
   input logic avm_local_bb4_ld__u24_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u24_inst0_writeack,
   // AVM avm_local_bb4_ld__u25_inst0
   output logic avm_local_bb4_ld__u25_inst0_enable,
   output logic avm_local_bb4_ld__u25_inst0_read,
   output logic avm_local_bb4_ld__u25_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u25_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u25_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u25_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u25_inst0_byteenable,
   input logic avm_local_bb4_ld__u25_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u25_inst0_readdata,
   input logic avm_local_bb4_ld__u25_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u25_inst0_writeack,
   // AVM avm_local_bb4_ld__u26_inst0
   output logic avm_local_bb4_ld__u26_inst0_enable,
   output logic avm_local_bb4_ld__u26_inst0_read,
   output logic avm_local_bb4_ld__u26_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u26_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u26_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u26_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u26_inst0_byteenable,
   input logic avm_local_bb4_ld__u26_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u26_inst0_readdata,
   input logic avm_local_bb4_ld__u26_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u26_inst0_writeack,
   // AVM avm_local_bb4_ld__u27_inst0
   output logic avm_local_bb4_ld__u27_inst0_enable,
   output logic avm_local_bb4_ld__u27_inst0_read,
   output logic avm_local_bb4_ld__u27_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u27_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u27_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u27_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u27_inst0_byteenable,
   input logic avm_local_bb4_ld__u27_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u27_inst0_readdata,
   input logic avm_local_bb4_ld__u27_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u27_inst0_writeack,
   // AVM avm_local_bb4_ld__u28_inst0
   output logic avm_local_bb4_ld__u28_inst0_enable,
   output logic avm_local_bb4_ld__u28_inst0_read,
   output logic avm_local_bb4_ld__u28_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u28_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u28_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u28_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u28_inst0_byteenable,
   input logic avm_local_bb4_ld__u28_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u28_inst0_readdata,
   input logic avm_local_bb4_ld__u28_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u28_inst0_writeack,
   // AVM avm_local_bb4_ld__u29_inst0
   output logic avm_local_bb4_ld__u29_inst0_enable,
   output logic avm_local_bb4_ld__u29_inst0_read,
   output logic avm_local_bb4_ld__u29_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u29_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u29_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u29_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u29_inst0_byteenable,
   input logic avm_local_bb4_ld__u29_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u29_inst0_readdata,
   input logic avm_local_bb4_ld__u29_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u29_inst0_writeack,
   // AVM avm_local_bb4_ld__u30_inst0
   output logic avm_local_bb4_ld__u30_inst0_enable,
   output logic avm_local_bb4_ld__u30_inst0_read,
   output logic avm_local_bb4_ld__u30_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u30_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u30_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u30_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u30_inst0_byteenable,
   input logic avm_local_bb4_ld__u30_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u30_inst0_readdata,
   input logic avm_local_bb4_ld__u30_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u30_inst0_writeack,
   // AVM avm_local_bb4_ld__u31_inst0
   output logic avm_local_bb4_ld__u31_inst0_enable,
   output logic avm_local_bb4_ld__u31_inst0_read,
   output logic avm_local_bb4_ld__u31_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u31_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u31_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u31_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u31_inst0_byteenable,
   input logic avm_local_bb4_ld__u31_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u31_inst0_readdata,
   input logic avm_local_bb4_ld__u31_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u31_inst0_writeack,
   // AVM avm_local_bb4_ld__u32_inst0
   output logic avm_local_bb4_ld__u32_inst0_enable,
   output logic avm_local_bb4_ld__u32_inst0_read,
   output logic avm_local_bb4_ld__u32_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u32_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u32_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u32_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u32_inst0_byteenable,
   input logic avm_local_bb4_ld__u32_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u32_inst0_readdata,
   input logic avm_local_bb4_ld__u32_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u32_inst0_writeack,
   // AVM avm_local_bb4_ld__u33_inst0
   output logic avm_local_bb4_ld__u33_inst0_enable,
   output logic avm_local_bb4_ld__u33_inst0_read,
   output logic avm_local_bb4_ld__u33_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u33_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u33_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u33_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u33_inst0_byteenable,
   input logic avm_local_bb4_ld__u33_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u33_inst0_readdata,
   input logic avm_local_bb4_ld__u33_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u33_inst0_writeack,
   // AVM avm_local_bb4_ld__u34_inst0
   output logic avm_local_bb4_ld__u34_inst0_enable,
   output logic avm_local_bb4_ld__u34_inst0_read,
   output logic avm_local_bb4_ld__u34_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u34_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u34_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u34_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u34_inst0_byteenable,
   input logic avm_local_bb4_ld__u34_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u34_inst0_readdata,
   input logic avm_local_bb4_ld__u34_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u34_inst0_writeack,
   // AVM avm_local_bb4_ld__u7_inst0
   output logic avm_local_bb4_ld__u7_inst0_enable,
   output logic avm_local_bb4_ld__u7_inst0_read,
   output logic avm_local_bb4_ld__u7_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u7_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u7_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u7_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u7_inst0_byteenable,
   input logic avm_local_bb4_ld__u7_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u7_inst0_readdata,
   input logic avm_local_bb4_ld__u7_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u7_inst0_writeack,
   // AVST avm_channel_id_chan_Conf2Intere_active_read
   input logic avm_channel_id_chan_Conf2Intere_active_read_valid,
   output logic avm_channel_id_chan_Conf2Intere_active_read_ready,
   input logic [7:0] avm_channel_id_chan_Conf2Intere_active_read_data,
   // AVST avm_channel_id_chan_Conf2Intere_cnt_read
   input logic avm_channel_id_chan_Conf2Intere_cnt_read_valid,
   output logic avm_channel_id_chan_Conf2Intere_cnt_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intere_cnt_read_data,
   // AVST avm_channel_id_chan_Conf2Intere_mode_read
   input logic avm_channel_id_chan_Conf2Intere_mode_read_valid,
   output logic avm_channel_id_chan_Conf2Intere_mode_read_ready,
   input logic [7:0] avm_channel_id_chan_Conf2Intere_mode_read_data,
   // AVST avm_channel_id_chan_Conf2Intere_x_read
   input logic avm_channel_id_chan_Conf2Intere_x_read_valid,
   output logic avm_channel_id_chan_Conf2Intere_x_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intere_x_read_data,
   // AVST avm_channel_id_chan_Conf2Intere_y_read
   input logic avm_channel_id_chan_Conf2Intere_y_read_valid,
   output logic avm_channel_id_chan_Conf2Intere_y_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intere_y_read_data,
   // AVST avm_channel_id_chan_Conf2Intere_z_read
   input logic avm_channel_id_chan_Conf2Intere_z_read_valid,
   output logic avm_channel_id_chan_Conf2Intere_z_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intere_z_read_data,
   // AVST avm_channel_id_chan_Intere2Store_active_write
   output logic avm_channel_id_chan_Intere2Store_active_write_valid,
   input logic avm_channel_id_chan_Intere2Store_active_write_ready,
   output logic [7:0] avm_channel_id_chan_Intere2Store_active_write_data,
   input logic avm_channel_id_chan_Intere2Store_active_write_almostfull,
   // AVST avm_channel_id_chan_Intere2Store_cnt_write
   output logic avm_channel_id_chan_Intere2Store_cnt_write_valid,
   input logic avm_channel_id_chan_Intere2Store_cnt_write_ready,
   output logic [31:0] avm_channel_id_chan_Intere2Store_cnt_write_data,
   input logic avm_channel_id_chan_Intere2Store_cnt_write_almostfull,
   // AVST avm_channel_id_chan_Intere2Store_intere_write
   output logic avm_channel_id_chan_Intere2Store_intere_write_valid,
   input logic avm_channel_id_chan_Intere2Store_intere_write_ready,
   output logic [31:0] avm_channel_id_chan_Intere2Store_intere_write_data,
   input logic avm_channel_id_chan_Intere2Store_intere_write_almostfull,
   // AVST avm_channel_id_chan_Intere2Store_mode_write
   output logic avm_channel_id_chan_Intere2Store_mode_write_valid,
   input logic avm_channel_id_chan_Intere2Store_mode_write_ready,
   output logic [7:0] avm_channel_id_chan_Intere2Store_mode_write_data,
   input logic avm_channel_id_chan_Intere2Store_mode_write_almostfull,
   // AVM p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid
);
   genvar __i;
   genvar __j;
   genvar __k;
   genvar __l;
   logic lmem_invalid_single_bit;
   logic [2:0] lmem_invalid_aspaces;
   logic local_avm_aspace11_enable [1][2];
   logic local_avm_aspace11_read [1][2];
   logic local_avm_aspace11_write [1][2];
   logic local_avm_aspace11_burstcount [1][2];
   logic [31:0] local_avm_aspace11_address [1][2];
   logic [31:0] local_avm_aspace11_writedata [1][2];
   logic [3:0] local_avm_aspace11_byteenable [1][2];
   logic local_avm_aspace11_waitrequest [1][2];
   logic [31:0] local_avm_aspace11_readdata [1][2];
   logic local_avm_aspace11_readdatavalid [1][2];
   logic local_avm_aspace11_writeack [1][2];
   logic local_avm_aspace13_enable [1][2];
   logic local_avm_aspace13_read [1][2];
   logic local_avm_aspace13_write [1][2];
   logic local_avm_aspace13_burstcount [1][2];
   logic [31:0] local_avm_aspace13_address [1][2];
   logic [31:0] local_avm_aspace13_writedata [1][2];
   logic [3:0] local_avm_aspace13_byteenable [1][2];
   logic local_avm_aspace13_waitrequest [1][2];
   logic [31:0] local_avm_aspace13_readdata [1][2];
   logic local_avm_aspace13_readdatavalid [1][2];
   logic local_avm_aspace13_writeack [1][2];
   logic local_avm_aspace15_enable [1][2];
   logic local_avm_aspace15_read [1][2];
   logic local_avm_aspace15_write [1][2];
   logic local_avm_aspace15_burstcount [1][2];
   logic [31:0] local_avm_aspace15_address [1][2];
   logic [31:0] local_avm_aspace15_writedata [1][2];
   logic [3:0] local_avm_aspace15_byteenable [1][2];
   logic local_avm_aspace15_waitrequest [1][2];
   logic [31:0] local_avm_aspace15_readdata [1][2];
   logic local_avm_aspace15_readdatavalid [1][2];
   logic local_avm_aspace15_writeack [1][2];

   // INST kernel of Krnl_InterE_function_wrapper
   Krnl_InterE_function_wrapper kernel
   (
      .local_router_hang(lmem_invalid_single_bit),
      .start(start),
      .kernel_arguments(kernel_arguments),
      .work_dim(work_dim),
      .global_offset_0(global_offset[0]),
      .global_offset_1(global_offset[1]),
      .global_offset_2(global_offset[2]),
      .kernel_valid_out(kernel_valid_out),
      .has_a_write_pending(has_a_write_pending),
      .has_a_lsu_active(has_a_lsu_active),
      .global_id_0(global_id[0]),
      .global_id_1(global_id[1]),
      .global_id_2(global_id[2]),
      .local_id_0(local_id[0]),
      .local_id_1(local_id[1]),
      .local_id_2(local_id[2]),
      .group_id_0(group_id[0]),
      .group_id_1(group_id[1]),
      .group_id_2(group_id[2]),
      .global_size_0(global_size[0]),
      .global_size_1(global_size[1]),
      .global_size_2(global_size[2]),
      .local_size_0(local_size[0]),
      .local_size_1(local_size[1]),
      .local_size_2(local_size[2]),
      .num_groups_0(num_groups[0]),
      .num_groups_1(num_groups[1]),
      .num_groups_2(num_groups[2]),
      .workgroup_size(workgroup_size),
      .kernel_stall_out(kernel_stall_out),
      .kernel_valid_in(kernel_valid_in),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack),
      // AVM avm_local_bb3_st__inst0
      .avm_local_bb3_st__inst0_enable(avm_local_bb3_st__inst0_enable),
      .avm_local_bb3_st__inst0_read(avm_local_bb3_st__inst0_read),
      .avm_local_bb3_st__inst0_write(avm_local_bb3_st__inst0_write),
      .avm_local_bb3_st__inst0_burstcount(avm_local_bb3_st__inst0_burstcount),
      .avm_local_bb3_st__inst0_address(avm_local_bb3_st__inst0_address),
      .avm_local_bb3_st__inst0_writedata(avm_local_bb3_st__inst0_writedata),
      .avm_local_bb3_st__inst0_byteenable(avm_local_bb3_st__inst0_byteenable),
      .avm_local_bb3_st__inst0_waitrequest(avm_local_bb3_st__inst0_waitrequest),
      .avm_local_bb3_st__inst0_readdata(avm_local_bb3_st__inst0_readdata),
      .avm_local_bb3_st__inst0_readdatavalid(avm_local_bb3_st__inst0_readdatavalid),
      .avm_local_bb3_st__inst0_writeack(avm_local_bb3_st__inst0_writeack),
      // AVM avm_local_bb4_ld__inst0
      .avm_local_bb4_ld__inst0_enable(avm_local_bb4_ld__inst0_enable),
      .avm_local_bb4_ld__inst0_read(avm_local_bb4_ld__inst0_read),
      .avm_local_bb4_ld__inst0_write(avm_local_bb4_ld__inst0_write),
      .avm_local_bb4_ld__inst0_burstcount(avm_local_bb4_ld__inst0_burstcount),
      .avm_local_bb4_ld__inst0_address(avm_local_bb4_ld__inst0_address),
      .avm_local_bb4_ld__inst0_writedata(avm_local_bb4_ld__inst0_writedata),
      .avm_local_bb4_ld__inst0_byteenable(avm_local_bb4_ld__inst0_byteenable),
      .avm_local_bb4_ld__inst0_waitrequest(avm_local_bb4_ld__inst0_waitrequest),
      .avm_local_bb4_ld__inst0_readdata(avm_local_bb4_ld__inst0_readdata),
      .avm_local_bb4_ld__inst0_readdatavalid(avm_local_bb4_ld__inst0_readdatavalid),
      .avm_local_bb4_ld__inst0_writeack(avm_local_bb4_ld__inst0_writeack),
      // AVM avm_local_bb4_ld__u11_inst0
      .avm_local_bb4_ld__u11_inst0_enable(avm_local_bb4_ld__u11_inst0_enable),
      .avm_local_bb4_ld__u11_inst0_read(avm_local_bb4_ld__u11_inst0_read),
      .avm_local_bb4_ld__u11_inst0_write(avm_local_bb4_ld__u11_inst0_write),
      .avm_local_bb4_ld__u11_inst0_burstcount(avm_local_bb4_ld__u11_inst0_burstcount),
      .avm_local_bb4_ld__u11_inst0_address(avm_local_bb4_ld__u11_inst0_address),
      .avm_local_bb4_ld__u11_inst0_writedata(avm_local_bb4_ld__u11_inst0_writedata),
      .avm_local_bb4_ld__u11_inst0_byteenable(avm_local_bb4_ld__u11_inst0_byteenable),
      .avm_local_bb4_ld__u11_inst0_waitrequest(avm_local_bb4_ld__u11_inst0_waitrequest),
      .avm_local_bb4_ld__u11_inst0_readdata(avm_local_bb4_ld__u11_inst0_readdata),
      .avm_local_bb4_ld__u11_inst0_readdatavalid(avm_local_bb4_ld__u11_inst0_readdatavalid),
      .avm_local_bb4_ld__u11_inst0_writeack(avm_local_bb4_ld__u11_inst0_writeack),
      // AVM avm_local_bb4_ld__u12_inst0
      .avm_local_bb4_ld__u12_inst0_enable(avm_local_bb4_ld__u12_inst0_enable),
      .avm_local_bb4_ld__u12_inst0_read(avm_local_bb4_ld__u12_inst0_read),
      .avm_local_bb4_ld__u12_inst0_write(avm_local_bb4_ld__u12_inst0_write),
      .avm_local_bb4_ld__u12_inst0_burstcount(avm_local_bb4_ld__u12_inst0_burstcount),
      .avm_local_bb4_ld__u12_inst0_address(avm_local_bb4_ld__u12_inst0_address),
      .avm_local_bb4_ld__u12_inst0_writedata(avm_local_bb4_ld__u12_inst0_writedata),
      .avm_local_bb4_ld__u12_inst0_byteenable(avm_local_bb4_ld__u12_inst0_byteenable),
      .avm_local_bb4_ld__u12_inst0_waitrequest(avm_local_bb4_ld__u12_inst0_waitrequest),
      .avm_local_bb4_ld__u12_inst0_readdata(avm_local_bb4_ld__u12_inst0_readdata),
      .avm_local_bb4_ld__u12_inst0_readdatavalid(avm_local_bb4_ld__u12_inst0_readdatavalid),
      .avm_local_bb4_ld__u12_inst0_writeack(avm_local_bb4_ld__u12_inst0_writeack),
      // AVM avm_local_bb4_ld__u13_inst0
      .avm_local_bb4_ld__u13_inst0_enable(avm_local_bb4_ld__u13_inst0_enable),
      .avm_local_bb4_ld__u13_inst0_read(avm_local_bb4_ld__u13_inst0_read),
      .avm_local_bb4_ld__u13_inst0_write(avm_local_bb4_ld__u13_inst0_write),
      .avm_local_bb4_ld__u13_inst0_burstcount(avm_local_bb4_ld__u13_inst0_burstcount),
      .avm_local_bb4_ld__u13_inst0_address(avm_local_bb4_ld__u13_inst0_address),
      .avm_local_bb4_ld__u13_inst0_writedata(avm_local_bb4_ld__u13_inst0_writedata),
      .avm_local_bb4_ld__u13_inst0_byteenable(avm_local_bb4_ld__u13_inst0_byteenable),
      .avm_local_bb4_ld__u13_inst0_waitrequest(avm_local_bb4_ld__u13_inst0_waitrequest),
      .avm_local_bb4_ld__u13_inst0_readdata(avm_local_bb4_ld__u13_inst0_readdata),
      .avm_local_bb4_ld__u13_inst0_readdatavalid(avm_local_bb4_ld__u13_inst0_readdatavalid),
      .avm_local_bb4_ld__u13_inst0_writeack(avm_local_bb4_ld__u13_inst0_writeack),
      // AVM avm_local_bb4_ld__u14_inst0
      .avm_local_bb4_ld__u14_inst0_enable(avm_local_bb4_ld__u14_inst0_enable),
      .avm_local_bb4_ld__u14_inst0_read(avm_local_bb4_ld__u14_inst0_read),
      .avm_local_bb4_ld__u14_inst0_write(avm_local_bb4_ld__u14_inst0_write),
      .avm_local_bb4_ld__u14_inst0_burstcount(avm_local_bb4_ld__u14_inst0_burstcount),
      .avm_local_bb4_ld__u14_inst0_address(avm_local_bb4_ld__u14_inst0_address),
      .avm_local_bb4_ld__u14_inst0_writedata(avm_local_bb4_ld__u14_inst0_writedata),
      .avm_local_bb4_ld__u14_inst0_byteenable(avm_local_bb4_ld__u14_inst0_byteenable),
      .avm_local_bb4_ld__u14_inst0_waitrequest(avm_local_bb4_ld__u14_inst0_waitrequest),
      .avm_local_bb4_ld__u14_inst0_readdata(avm_local_bb4_ld__u14_inst0_readdata),
      .avm_local_bb4_ld__u14_inst0_readdatavalid(avm_local_bb4_ld__u14_inst0_readdatavalid),
      .avm_local_bb4_ld__u14_inst0_writeack(avm_local_bb4_ld__u14_inst0_writeack),
      // AVM avm_local_bb4_ld__u15_inst0
      .avm_local_bb4_ld__u15_inst0_enable(avm_local_bb4_ld__u15_inst0_enable),
      .avm_local_bb4_ld__u15_inst0_read(avm_local_bb4_ld__u15_inst0_read),
      .avm_local_bb4_ld__u15_inst0_write(avm_local_bb4_ld__u15_inst0_write),
      .avm_local_bb4_ld__u15_inst0_burstcount(avm_local_bb4_ld__u15_inst0_burstcount),
      .avm_local_bb4_ld__u15_inst0_address(avm_local_bb4_ld__u15_inst0_address),
      .avm_local_bb4_ld__u15_inst0_writedata(avm_local_bb4_ld__u15_inst0_writedata),
      .avm_local_bb4_ld__u15_inst0_byteenable(avm_local_bb4_ld__u15_inst0_byteenable),
      .avm_local_bb4_ld__u15_inst0_waitrequest(avm_local_bb4_ld__u15_inst0_waitrequest),
      .avm_local_bb4_ld__u15_inst0_readdata(avm_local_bb4_ld__u15_inst0_readdata),
      .avm_local_bb4_ld__u15_inst0_readdatavalid(avm_local_bb4_ld__u15_inst0_readdatavalid),
      .avm_local_bb4_ld__u15_inst0_writeack(avm_local_bb4_ld__u15_inst0_writeack),
      // AVM avm_local_bb4_ld__u16_inst0
      .avm_local_bb4_ld__u16_inst0_enable(avm_local_bb4_ld__u16_inst0_enable),
      .avm_local_bb4_ld__u16_inst0_read(avm_local_bb4_ld__u16_inst0_read),
      .avm_local_bb4_ld__u16_inst0_write(avm_local_bb4_ld__u16_inst0_write),
      .avm_local_bb4_ld__u16_inst0_burstcount(avm_local_bb4_ld__u16_inst0_burstcount),
      .avm_local_bb4_ld__u16_inst0_address(avm_local_bb4_ld__u16_inst0_address),
      .avm_local_bb4_ld__u16_inst0_writedata(avm_local_bb4_ld__u16_inst0_writedata),
      .avm_local_bb4_ld__u16_inst0_byteenable(avm_local_bb4_ld__u16_inst0_byteenable),
      .avm_local_bb4_ld__u16_inst0_waitrequest(avm_local_bb4_ld__u16_inst0_waitrequest),
      .avm_local_bb4_ld__u16_inst0_readdata(avm_local_bb4_ld__u16_inst0_readdata),
      .avm_local_bb4_ld__u16_inst0_readdatavalid(avm_local_bb4_ld__u16_inst0_readdatavalid),
      .avm_local_bb4_ld__u16_inst0_writeack(avm_local_bb4_ld__u16_inst0_writeack),
      // AVM avm_local_bb4_ld__u17_inst0
      .avm_local_bb4_ld__u17_inst0_enable(avm_local_bb4_ld__u17_inst0_enable),
      .avm_local_bb4_ld__u17_inst0_read(avm_local_bb4_ld__u17_inst0_read),
      .avm_local_bb4_ld__u17_inst0_write(avm_local_bb4_ld__u17_inst0_write),
      .avm_local_bb4_ld__u17_inst0_burstcount(avm_local_bb4_ld__u17_inst0_burstcount),
      .avm_local_bb4_ld__u17_inst0_address(avm_local_bb4_ld__u17_inst0_address),
      .avm_local_bb4_ld__u17_inst0_writedata(avm_local_bb4_ld__u17_inst0_writedata),
      .avm_local_bb4_ld__u17_inst0_byteenable(avm_local_bb4_ld__u17_inst0_byteenable),
      .avm_local_bb4_ld__u17_inst0_waitrequest(avm_local_bb4_ld__u17_inst0_waitrequest),
      .avm_local_bb4_ld__u17_inst0_readdata(avm_local_bb4_ld__u17_inst0_readdata),
      .avm_local_bb4_ld__u17_inst0_readdatavalid(avm_local_bb4_ld__u17_inst0_readdatavalid),
      .avm_local_bb4_ld__u17_inst0_writeack(avm_local_bb4_ld__u17_inst0_writeack),
      // AVM avm_local_bb4_ld__u18_inst0
      .avm_local_bb4_ld__u18_inst0_enable(avm_local_bb4_ld__u18_inst0_enable),
      .avm_local_bb4_ld__u18_inst0_read(avm_local_bb4_ld__u18_inst0_read),
      .avm_local_bb4_ld__u18_inst0_write(avm_local_bb4_ld__u18_inst0_write),
      .avm_local_bb4_ld__u18_inst0_burstcount(avm_local_bb4_ld__u18_inst0_burstcount),
      .avm_local_bb4_ld__u18_inst0_address(avm_local_bb4_ld__u18_inst0_address),
      .avm_local_bb4_ld__u18_inst0_writedata(avm_local_bb4_ld__u18_inst0_writedata),
      .avm_local_bb4_ld__u18_inst0_byteenable(avm_local_bb4_ld__u18_inst0_byteenable),
      .avm_local_bb4_ld__u18_inst0_waitrequest(avm_local_bb4_ld__u18_inst0_waitrequest),
      .avm_local_bb4_ld__u18_inst0_readdata(avm_local_bb4_ld__u18_inst0_readdata),
      .avm_local_bb4_ld__u18_inst0_readdatavalid(avm_local_bb4_ld__u18_inst0_readdatavalid),
      .avm_local_bb4_ld__u18_inst0_writeack(avm_local_bb4_ld__u18_inst0_writeack),
      // AVM avm_local_bb4_ld__u19_inst0
      .avm_local_bb4_ld__u19_inst0_enable(avm_local_bb4_ld__u19_inst0_enable),
      .avm_local_bb4_ld__u19_inst0_read(avm_local_bb4_ld__u19_inst0_read),
      .avm_local_bb4_ld__u19_inst0_write(avm_local_bb4_ld__u19_inst0_write),
      .avm_local_bb4_ld__u19_inst0_burstcount(avm_local_bb4_ld__u19_inst0_burstcount),
      .avm_local_bb4_ld__u19_inst0_address(avm_local_bb4_ld__u19_inst0_address),
      .avm_local_bb4_ld__u19_inst0_writedata(avm_local_bb4_ld__u19_inst0_writedata),
      .avm_local_bb4_ld__u19_inst0_byteenable(avm_local_bb4_ld__u19_inst0_byteenable),
      .avm_local_bb4_ld__u19_inst0_waitrequest(avm_local_bb4_ld__u19_inst0_waitrequest),
      .avm_local_bb4_ld__u19_inst0_readdata(avm_local_bb4_ld__u19_inst0_readdata),
      .avm_local_bb4_ld__u19_inst0_readdatavalid(avm_local_bb4_ld__u19_inst0_readdatavalid),
      .avm_local_bb4_ld__u19_inst0_writeack(avm_local_bb4_ld__u19_inst0_writeack),
      // AVM avm_local_bb4_ld__u20_inst0
      .avm_local_bb4_ld__u20_inst0_enable(avm_local_bb4_ld__u20_inst0_enable),
      .avm_local_bb4_ld__u20_inst0_read(avm_local_bb4_ld__u20_inst0_read),
      .avm_local_bb4_ld__u20_inst0_write(avm_local_bb4_ld__u20_inst0_write),
      .avm_local_bb4_ld__u20_inst0_burstcount(avm_local_bb4_ld__u20_inst0_burstcount),
      .avm_local_bb4_ld__u20_inst0_address(avm_local_bb4_ld__u20_inst0_address),
      .avm_local_bb4_ld__u20_inst0_writedata(avm_local_bb4_ld__u20_inst0_writedata),
      .avm_local_bb4_ld__u20_inst0_byteenable(avm_local_bb4_ld__u20_inst0_byteenable),
      .avm_local_bb4_ld__u20_inst0_waitrequest(avm_local_bb4_ld__u20_inst0_waitrequest),
      .avm_local_bb4_ld__u20_inst0_readdata(avm_local_bb4_ld__u20_inst0_readdata),
      .avm_local_bb4_ld__u20_inst0_readdatavalid(avm_local_bb4_ld__u20_inst0_readdatavalid),
      .avm_local_bb4_ld__u20_inst0_writeack(avm_local_bb4_ld__u20_inst0_writeack),
      // AVM avm_local_bb4_ld__u21_inst0
      .avm_local_bb4_ld__u21_inst0_enable(avm_local_bb4_ld__u21_inst0_enable),
      .avm_local_bb4_ld__u21_inst0_read(avm_local_bb4_ld__u21_inst0_read),
      .avm_local_bb4_ld__u21_inst0_write(avm_local_bb4_ld__u21_inst0_write),
      .avm_local_bb4_ld__u21_inst0_burstcount(avm_local_bb4_ld__u21_inst0_burstcount),
      .avm_local_bb4_ld__u21_inst0_address(avm_local_bb4_ld__u21_inst0_address),
      .avm_local_bb4_ld__u21_inst0_writedata(avm_local_bb4_ld__u21_inst0_writedata),
      .avm_local_bb4_ld__u21_inst0_byteenable(avm_local_bb4_ld__u21_inst0_byteenable),
      .avm_local_bb4_ld__u21_inst0_waitrequest(avm_local_bb4_ld__u21_inst0_waitrequest),
      .avm_local_bb4_ld__u21_inst0_readdata(avm_local_bb4_ld__u21_inst0_readdata),
      .avm_local_bb4_ld__u21_inst0_readdatavalid(avm_local_bb4_ld__u21_inst0_readdatavalid),
      .avm_local_bb4_ld__u21_inst0_writeack(avm_local_bb4_ld__u21_inst0_writeack),
      // AVM avm_local_bb4_ld__u22_inst0
      .avm_local_bb4_ld__u22_inst0_enable(avm_local_bb4_ld__u22_inst0_enable),
      .avm_local_bb4_ld__u22_inst0_read(avm_local_bb4_ld__u22_inst0_read),
      .avm_local_bb4_ld__u22_inst0_write(avm_local_bb4_ld__u22_inst0_write),
      .avm_local_bb4_ld__u22_inst0_burstcount(avm_local_bb4_ld__u22_inst0_burstcount),
      .avm_local_bb4_ld__u22_inst0_address(avm_local_bb4_ld__u22_inst0_address),
      .avm_local_bb4_ld__u22_inst0_writedata(avm_local_bb4_ld__u22_inst0_writedata),
      .avm_local_bb4_ld__u22_inst0_byteenable(avm_local_bb4_ld__u22_inst0_byteenable),
      .avm_local_bb4_ld__u22_inst0_waitrequest(avm_local_bb4_ld__u22_inst0_waitrequest),
      .avm_local_bb4_ld__u22_inst0_readdata(avm_local_bb4_ld__u22_inst0_readdata),
      .avm_local_bb4_ld__u22_inst0_readdatavalid(avm_local_bb4_ld__u22_inst0_readdatavalid),
      .avm_local_bb4_ld__u22_inst0_writeack(avm_local_bb4_ld__u22_inst0_writeack),
      // AVM avm_local_bb4_ld__u23_inst0
      .avm_local_bb4_ld__u23_inst0_enable(avm_local_bb4_ld__u23_inst0_enable),
      .avm_local_bb4_ld__u23_inst0_read(avm_local_bb4_ld__u23_inst0_read),
      .avm_local_bb4_ld__u23_inst0_write(avm_local_bb4_ld__u23_inst0_write),
      .avm_local_bb4_ld__u23_inst0_burstcount(avm_local_bb4_ld__u23_inst0_burstcount),
      .avm_local_bb4_ld__u23_inst0_address(avm_local_bb4_ld__u23_inst0_address),
      .avm_local_bb4_ld__u23_inst0_writedata(avm_local_bb4_ld__u23_inst0_writedata),
      .avm_local_bb4_ld__u23_inst0_byteenable(avm_local_bb4_ld__u23_inst0_byteenable),
      .avm_local_bb4_ld__u23_inst0_waitrequest(avm_local_bb4_ld__u23_inst0_waitrequest),
      .avm_local_bb4_ld__u23_inst0_readdata(avm_local_bb4_ld__u23_inst0_readdata),
      .avm_local_bb4_ld__u23_inst0_readdatavalid(avm_local_bb4_ld__u23_inst0_readdatavalid),
      .avm_local_bb4_ld__u23_inst0_writeack(avm_local_bb4_ld__u23_inst0_writeack),
      // AVM avm_local_bb4_ld__u24_inst0
      .avm_local_bb4_ld__u24_inst0_enable(avm_local_bb4_ld__u24_inst0_enable),
      .avm_local_bb4_ld__u24_inst0_read(avm_local_bb4_ld__u24_inst0_read),
      .avm_local_bb4_ld__u24_inst0_write(avm_local_bb4_ld__u24_inst0_write),
      .avm_local_bb4_ld__u24_inst0_burstcount(avm_local_bb4_ld__u24_inst0_burstcount),
      .avm_local_bb4_ld__u24_inst0_address(avm_local_bb4_ld__u24_inst0_address),
      .avm_local_bb4_ld__u24_inst0_writedata(avm_local_bb4_ld__u24_inst0_writedata),
      .avm_local_bb4_ld__u24_inst0_byteenable(avm_local_bb4_ld__u24_inst0_byteenable),
      .avm_local_bb4_ld__u24_inst0_waitrequest(avm_local_bb4_ld__u24_inst0_waitrequest),
      .avm_local_bb4_ld__u24_inst0_readdata(avm_local_bb4_ld__u24_inst0_readdata),
      .avm_local_bb4_ld__u24_inst0_readdatavalid(avm_local_bb4_ld__u24_inst0_readdatavalid),
      .avm_local_bb4_ld__u24_inst0_writeack(avm_local_bb4_ld__u24_inst0_writeack),
      // AVM avm_local_bb4_ld__u25_inst0
      .avm_local_bb4_ld__u25_inst0_enable(avm_local_bb4_ld__u25_inst0_enable),
      .avm_local_bb4_ld__u25_inst0_read(avm_local_bb4_ld__u25_inst0_read),
      .avm_local_bb4_ld__u25_inst0_write(avm_local_bb4_ld__u25_inst0_write),
      .avm_local_bb4_ld__u25_inst0_burstcount(avm_local_bb4_ld__u25_inst0_burstcount),
      .avm_local_bb4_ld__u25_inst0_address(avm_local_bb4_ld__u25_inst0_address),
      .avm_local_bb4_ld__u25_inst0_writedata(avm_local_bb4_ld__u25_inst0_writedata),
      .avm_local_bb4_ld__u25_inst0_byteenable(avm_local_bb4_ld__u25_inst0_byteenable),
      .avm_local_bb4_ld__u25_inst0_waitrequest(avm_local_bb4_ld__u25_inst0_waitrequest),
      .avm_local_bb4_ld__u25_inst0_readdata(avm_local_bb4_ld__u25_inst0_readdata),
      .avm_local_bb4_ld__u25_inst0_readdatavalid(avm_local_bb4_ld__u25_inst0_readdatavalid),
      .avm_local_bb4_ld__u25_inst0_writeack(avm_local_bb4_ld__u25_inst0_writeack),
      // AVM avm_local_bb4_ld__u26_inst0
      .avm_local_bb4_ld__u26_inst0_enable(avm_local_bb4_ld__u26_inst0_enable),
      .avm_local_bb4_ld__u26_inst0_read(avm_local_bb4_ld__u26_inst0_read),
      .avm_local_bb4_ld__u26_inst0_write(avm_local_bb4_ld__u26_inst0_write),
      .avm_local_bb4_ld__u26_inst0_burstcount(avm_local_bb4_ld__u26_inst0_burstcount),
      .avm_local_bb4_ld__u26_inst0_address(avm_local_bb4_ld__u26_inst0_address),
      .avm_local_bb4_ld__u26_inst0_writedata(avm_local_bb4_ld__u26_inst0_writedata),
      .avm_local_bb4_ld__u26_inst0_byteenable(avm_local_bb4_ld__u26_inst0_byteenable),
      .avm_local_bb4_ld__u26_inst0_waitrequest(avm_local_bb4_ld__u26_inst0_waitrequest),
      .avm_local_bb4_ld__u26_inst0_readdata(avm_local_bb4_ld__u26_inst0_readdata),
      .avm_local_bb4_ld__u26_inst0_readdatavalid(avm_local_bb4_ld__u26_inst0_readdatavalid),
      .avm_local_bb4_ld__u26_inst0_writeack(avm_local_bb4_ld__u26_inst0_writeack),
      // AVM avm_local_bb4_ld__u27_inst0
      .avm_local_bb4_ld__u27_inst0_enable(avm_local_bb4_ld__u27_inst0_enable),
      .avm_local_bb4_ld__u27_inst0_read(avm_local_bb4_ld__u27_inst0_read),
      .avm_local_bb4_ld__u27_inst0_write(avm_local_bb4_ld__u27_inst0_write),
      .avm_local_bb4_ld__u27_inst0_burstcount(avm_local_bb4_ld__u27_inst0_burstcount),
      .avm_local_bb4_ld__u27_inst0_address(avm_local_bb4_ld__u27_inst0_address),
      .avm_local_bb4_ld__u27_inst0_writedata(avm_local_bb4_ld__u27_inst0_writedata),
      .avm_local_bb4_ld__u27_inst0_byteenable(avm_local_bb4_ld__u27_inst0_byteenable),
      .avm_local_bb4_ld__u27_inst0_waitrequest(avm_local_bb4_ld__u27_inst0_waitrequest),
      .avm_local_bb4_ld__u27_inst0_readdata(avm_local_bb4_ld__u27_inst0_readdata),
      .avm_local_bb4_ld__u27_inst0_readdatavalid(avm_local_bb4_ld__u27_inst0_readdatavalid),
      .avm_local_bb4_ld__u27_inst0_writeack(avm_local_bb4_ld__u27_inst0_writeack),
      // AVM avm_local_bb4_ld__u28_inst0
      .avm_local_bb4_ld__u28_inst0_enable(avm_local_bb4_ld__u28_inst0_enable),
      .avm_local_bb4_ld__u28_inst0_read(avm_local_bb4_ld__u28_inst0_read),
      .avm_local_bb4_ld__u28_inst0_write(avm_local_bb4_ld__u28_inst0_write),
      .avm_local_bb4_ld__u28_inst0_burstcount(avm_local_bb4_ld__u28_inst0_burstcount),
      .avm_local_bb4_ld__u28_inst0_address(avm_local_bb4_ld__u28_inst0_address),
      .avm_local_bb4_ld__u28_inst0_writedata(avm_local_bb4_ld__u28_inst0_writedata),
      .avm_local_bb4_ld__u28_inst0_byteenable(avm_local_bb4_ld__u28_inst0_byteenable),
      .avm_local_bb4_ld__u28_inst0_waitrequest(avm_local_bb4_ld__u28_inst0_waitrequest),
      .avm_local_bb4_ld__u28_inst0_readdata(avm_local_bb4_ld__u28_inst0_readdata),
      .avm_local_bb4_ld__u28_inst0_readdatavalid(avm_local_bb4_ld__u28_inst0_readdatavalid),
      .avm_local_bb4_ld__u28_inst0_writeack(avm_local_bb4_ld__u28_inst0_writeack),
      // AVM avm_local_bb4_ld__u29_inst0
      .avm_local_bb4_ld__u29_inst0_enable(avm_local_bb4_ld__u29_inst0_enable),
      .avm_local_bb4_ld__u29_inst0_read(avm_local_bb4_ld__u29_inst0_read),
      .avm_local_bb4_ld__u29_inst0_write(avm_local_bb4_ld__u29_inst0_write),
      .avm_local_bb4_ld__u29_inst0_burstcount(avm_local_bb4_ld__u29_inst0_burstcount),
      .avm_local_bb4_ld__u29_inst0_address(avm_local_bb4_ld__u29_inst0_address),
      .avm_local_bb4_ld__u29_inst0_writedata(avm_local_bb4_ld__u29_inst0_writedata),
      .avm_local_bb4_ld__u29_inst0_byteenable(avm_local_bb4_ld__u29_inst0_byteenable),
      .avm_local_bb4_ld__u29_inst0_waitrequest(avm_local_bb4_ld__u29_inst0_waitrequest),
      .avm_local_bb4_ld__u29_inst0_readdata(avm_local_bb4_ld__u29_inst0_readdata),
      .avm_local_bb4_ld__u29_inst0_readdatavalid(avm_local_bb4_ld__u29_inst0_readdatavalid),
      .avm_local_bb4_ld__u29_inst0_writeack(avm_local_bb4_ld__u29_inst0_writeack),
      // AVM avm_local_bb4_ld__u30_inst0
      .avm_local_bb4_ld__u30_inst0_enable(avm_local_bb4_ld__u30_inst0_enable),
      .avm_local_bb4_ld__u30_inst0_read(avm_local_bb4_ld__u30_inst0_read),
      .avm_local_bb4_ld__u30_inst0_write(avm_local_bb4_ld__u30_inst0_write),
      .avm_local_bb4_ld__u30_inst0_burstcount(avm_local_bb4_ld__u30_inst0_burstcount),
      .avm_local_bb4_ld__u30_inst0_address(avm_local_bb4_ld__u30_inst0_address),
      .avm_local_bb4_ld__u30_inst0_writedata(avm_local_bb4_ld__u30_inst0_writedata),
      .avm_local_bb4_ld__u30_inst0_byteenable(avm_local_bb4_ld__u30_inst0_byteenable),
      .avm_local_bb4_ld__u30_inst0_waitrequest(avm_local_bb4_ld__u30_inst0_waitrequest),
      .avm_local_bb4_ld__u30_inst0_readdata(avm_local_bb4_ld__u30_inst0_readdata),
      .avm_local_bb4_ld__u30_inst0_readdatavalid(avm_local_bb4_ld__u30_inst0_readdatavalid),
      .avm_local_bb4_ld__u30_inst0_writeack(avm_local_bb4_ld__u30_inst0_writeack),
      // AVM avm_local_bb4_ld__u31_inst0
      .avm_local_bb4_ld__u31_inst0_enable(avm_local_bb4_ld__u31_inst0_enable),
      .avm_local_bb4_ld__u31_inst0_read(avm_local_bb4_ld__u31_inst0_read),
      .avm_local_bb4_ld__u31_inst0_write(avm_local_bb4_ld__u31_inst0_write),
      .avm_local_bb4_ld__u31_inst0_burstcount(avm_local_bb4_ld__u31_inst0_burstcount),
      .avm_local_bb4_ld__u31_inst0_address(avm_local_bb4_ld__u31_inst0_address),
      .avm_local_bb4_ld__u31_inst0_writedata(avm_local_bb4_ld__u31_inst0_writedata),
      .avm_local_bb4_ld__u31_inst0_byteenable(avm_local_bb4_ld__u31_inst0_byteenable),
      .avm_local_bb4_ld__u31_inst0_waitrequest(avm_local_bb4_ld__u31_inst0_waitrequest),
      .avm_local_bb4_ld__u31_inst0_readdata(avm_local_bb4_ld__u31_inst0_readdata),
      .avm_local_bb4_ld__u31_inst0_readdatavalid(avm_local_bb4_ld__u31_inst0_readdatavalid),
      .avm_local_bb4_ld__u31_inst0_writeack(avm_local_bb4_ld__u31_inst0_writeack),
      // AVM avm_local_bb4_ld__u32_inst0
      .avm_local_bb4_ld__u32_inst0_enable(avm_local_bb4_ld__u32_inst0_enable),
      .avm_local_bb4_ld__u32_inst0_read(avm_local_bb4_ld__u32_inst0_read),
      .avm_local_bb4_ld__u32_inst0_write(avm_local_bb4_ld__u32_inst0_write),
      .avm_local_bb4_ld__u32_inst0_burstcount(avm_local_bb4_ld__u32_inst0_burstcount),
      .avm_local_bb4_ld__u32_inst0_address(avm_local_bb4_ld__u32_inst0_address),
      .avm_local_bb4_ld__u32_inst0_writedata(avm_local_bb4_ld__u32_inst0_writedata),
      .avm_local_bb4_ld__u32_inst0_byteenable(avm_local_bb4_ld__u32_inst0_byteenable),
      .avm_local_bb4_ld__u32_inst0_waitrequest(avm_local_bb4_ld__u32_inst0_waitrequest),
      .avm_local_bb4_ld__u32_inst0_readdata(avm_local_bb4_ld__u32_inst0_readdata),
      .avm_local_bb4_ld__u32_inst0_readdatavalid(avm_local_bb4_ld__u32_inst0_readdatavalid),
      .avm_local_bb4_ld__u32_inst0_writeack(avm_local_bb4_ld__u32_inst0_writeack),
      // AVM avm_local_bb4_ld__u33_inst0
      .avm_local_bb4_ld__u33_inst0_enable(avm_local_bb4_ld__u33_inst0_enable),
      .avm_local_bb4_ld__u33_inst0_read(avm_local_bb4_ld__u33_inst0_read),
      .avm_local_bb4_ld__u33_inst0_write(avm_local_bb4_ld__u33_inst0_write),
      .avm_local_bb4_ld__u33_inst0_burstcount(avm_local_bb4_ld__u33_inst0_burstcount),
      .avm_local_bb4_ld__u33_inst0_address(avm_local_bb4_ld__u33_inst0_address),
      .avm_local_bb4_ld__u33_inst0_writedata(avm_local_bb4_ld__u33_inst0_writedata),
      .avm_local_bb4_ld__u33_inst0_byteenable(avm_local_bb4_ld__u33_inst0_byteenable),
      .avm_local_bb4_ld__u33_inst0_waitrequest(avm_local_bb4_ld__u33_inst0_waitrequest),
      .avm_local_bb4_ld__u33_inst0_readdata(avm_local_bb4_ld__u33_inst0_readdata),
      .avm_local_bb4_ld__u33_inst0_readdatavalid(avm_local_bb4_ld__u33_inst0_readdatavalid),
      .avm_local_bb4_ld__u33_inst0_writeack(avm_local_bb4_ld__u33_inst0_writeack),
      // AVM avm_local_bb4_ld__u34_inst0
      .avm_local_bb4_ld__u34_inst0_enable(avm_local_bb4_ld__u34_inst0_enable),
      .avm_local_bb4_ld__u34_inst0_read(avm_local_bb4_ld__u34_inst0_read),
      .avm_local_bb4_ld__u34_inst0_write(avm_local_bb4_ld__u34_inst0_write),
      .avm_local_bb4_ld__u34_inst0_burstcount(avm_local_bb4_ld__u34_inst0_burstcount),
      .avm_local_bb4_ld__u34_inst0_address(avm_local_bb4_ld__u34_inst0_address),
      .avm_local_bb4_ld__u34_inst0_writedata(avm_local_bb4_ld__u34_inst0_writedata),
      .avm_local_bb4_ld__u34_inst0_byteenable(avm_local_bb4_ld__u34_inst0_byteenable),
      .avm_local_bb4_ld__u34_inst0_waitrequest(avm_local_bb4_ld__u34_inst0_waitrequest),
      .avm_local_bb4_ld__u34_inst0_readdata(avm_local_bb4_ld__u34_inst0_readdata),
      .avm_local_bb4_ld__u34_inst0_readdatavalid(avm_local_bb4_ld__u34_inst0_readdatavalid),
      .avm_local_bb4_ld__u34_inst0_writeack(avm_local_bb4_ld__u34_inst0_writeack),
      // AVM avm_local_bb4_ld__u7_inst0
      .avm_local_bb4_ld__u7_inst0_enable(avm_local_bb4_ld__u7_inst0_enable),
      .avm_local_bb4_ld__u7_inst0_read(avm_local_bb4_ld__u7_inst0_read),
      .avm_local_bb4_ld__u7_inst0_write(avm_local_bb4_ld__u7_inst0_write),
      .avm_local_bb4_ld__u7_inst0_burstcount(avm_local_bb4_ld__u7_inst0_burstcount),
      .avm_local_bb4_ld__u7_inst0_address(avm_local_bb4_ld__u7_inst0_address),
      .avm_local_bb4_ld__u7_inst0_writedata(avm_local_bb4_ld__u7_inst0_writedata),
      .avm_local_bb4_ld__u7_inst0_byteenable(avm_local_bb4_ld__u7_inst0_byteenable),
      .avm_local_bb4_ld__u7_inst0_waitrequest(avm_local_bb4_ld__u7_inst0_waitrequest),
      .avm_local_bb4_ld__u7_inst0_readdata(avm_local_bb4_ld__u7_inst0_readdata),
      .avm_local_bb4_ld__u7_inst0_readdatavalid(avm_local_bb4_ld__u7_inst0_readdatavalid),
      .avm_local_bb4_ld__u7_inst0_writeack(avm_local_bb4_ld__u7_inst0_writeack),
      // AVM avm_local_bb2_st__inst0
      .avm_local_bb2_st__inst0_enable(local_avm_aspace11_enable[0][0]),
      .avm_local_bb2_st__inst0_read(local_avm_aspace11_read[0][0]),
      .avm_local_bb2_st__inst0_write(local_avm_aspace11_write[0][0]),
      .avm_local_bb2_st__inst0_burstcount(local_avm_aspace11_burstcount[0][0]),
      .avm_local_bb2_st__inst0_address(local_avm_aspace11_address[0][0]),
      .avm_local_bb2_st__inst0_writedata(local_avm_aspace11_writedata[0][0]),
      .avm_local_bb2_st__inst0_byteenable(local_avm_aspace11_byteenable[0][0]),
      .avm_local_bb2_st__inst0_waitrequest(local_avm_aspace11_waitrequest[0][0]),
      .avm_local_bb2_st__inst0_readdata(local_avm_aspace11_readdata[0][0]),
      .avm_local_bb2_st__inst0_readdatavalid(local_avm_aspace11_readdatavalid[0][0]),
      .avm_local_bb2_st__inst0_writeack(local_avm_aspace11_writeack[0][0]),
      // AVM avm_local_bb4_ld__u9_inst0
      .avm_local_bb4_ld__u9_inst0_enable(local_avm_aspace11_enable[0][1]),
      .avm_local_bb4_ld__u9_inst0_read(local_avm_aspace11_read[0][1]),
      .avm_local_bb4_ld__u9_inst0_write(local_avm_aspace11_write[0][1]),
      .avm_local_bb4_ld__u9_inst0_burstcount(local_avm_aspace11_burstcount[0][1]),
      .avm_local_bb4_ld__u9_inst0_address(local_avm_aspace11_address[0][1]),
      .avm_local_bb4_ld__u9_inst0_writedata(local_avm_aspace11_writedata[0][1]),
      .avm_local_bb4_ld__u9_inst0_byteenable(local_avm_aspace11_byteenable[0][1]),
      .avm_local_bb4_ld__u9_inst0_waitrequest(local_avm_aspace11_waitrequest[0][1]),
      .avm_local_bb4_ld__u9_inst0_readdata(local_avm_aspace11_readdata[0][1]),
      .avm_local_bb4_ld__u9_inst0_readdatavalid(local_avm_aspace11_readdatavalid[0][1]),
      .avm_local_bb4_ld__u9_inst0_writeack(local_avm_aspace11_writeack[0][1]),
      // AVM avm_local_bb2_st__u0_inst0
      .avm_local_bb2_st__u0_inst0_enable(local_avm_aspace13_enable[0][0]),
      .avm_local_bb2_st__u0_inst0_read(local_avm_aspace13_read[0][0]),
      .avm_local_bb2_st__u0_inst0_write(local_avm_aspace13_write[0][0]),
      .avm_local_bb2_st__u0_inst0_burstcount(local_avm_aspace13_burstcount[0][0]),
      .avm_local_bb2_st__u0_inst0_address(local_avm_aspace13_address[0][0]),
      .avm_local_bb2_st__u0_inst0_writedata(local_avm_aspace13_writedata[0][0]),
      .avm_local_bb2_st__u0_inst0_byteenable(local_avm_aspace13_byteenable[0][0]),
      .avm_local_bb2_st__u0_inst0_waitrequest(local_avm_aspace13_waitrequest[0][0]),
      .avm_local_bb2_st__u0_inst0_readdata(local_avm_aspace13_readdata[0][0]),
      .avm_local_bb2_st__u0_inst0_readdatavalid(local_avm_aspace13_readdatavalid[0][0]),
      .avm_local_bb2_st__u0_inst0_writeack(local_avm_aspace13_writeack[0][0]),
      // AVM avm_local_bb4_ld__u8_inst0
      .avm_local_bb4_ld__u8_inst0_enable(local_avm_aspace13_enable[0][1]),
      .avm_local_bb4_ld__u8_inst0_read(local_avm_aspace13_read[0][1]),
      .avm_local_bb4_ld__u8_inst0_write(local_avm_aspace13_write[0][1]),
      .avm_local_bb4_ld__u8_inst0_burstcount(local_avm_aspace13_burstcount[0][1]),
      .avm_local_bb4_ld__u8_inst0_address(local_avm_aspace13_address[0][1]),
      .avm_local_bb4_ld__u8_inst0_writedata(local_avm_aspace13_writedata[0][1]),
      .avm_local_bb4_ld__u8_inst0_byteenable(local_avm_aspace13_byteenable[0][1]),
      .avm_local_bb4_ld__u8_inst0_waitrequest(local_avm_aspace13_waitrequest[0][1]),
      .avm_local_bb4_ld__u8_inst0_readdata(local_avm_aspace13_readdata[0][1]),
      .avm_local_bb4_ld__u8_inst0_readdatavalid(local_avm_aspace13_readdatavalid[0][1]),
      .avm_local_bb4_ld__u8_inst0_writeack(local_avm_aspace13_writeack[0][1]),
      // AVM avm_local_bb2_st__u2_inst0
      .avm_local_bb2_st__u2_inst0_enable(local_avm_aspace15_enable[0][0]),
      .avm_local_bb2_st__u2_inst0_read(local_avm_aspace15_read[0][0]),
      .avm_local_bb2_st__u2_inst0_write(local_avm_aspace15_write[0][0]),
      .avm_local_bb2_st__u2_inst0_burstcount(local_avm_aspace15_burstcount[0][0]),
      .avm_local_bb2_st__u2_inst0_address(local_avm_aspace15_address[0][0]),
      .avm_local_bb2_st__u2_inst0_writedata(local_avm_aspace15_writedata[0][0]),
      .avm_local_bb2_st__u2_inst0_byteenable(local_avm_aspace15_byteenable[0][0]),
      .avm_local_bb2_st__u2_inst0_waitrequest(local_avm_aspace15_waitrequest[0][0]),
      .avm_local_bb2_st__u2_inst0_readdata(local_avm_aspace15_readdata[0][0]),
      .avm_local_bb2_st__u2_inst0_readdatavalid(local_avm_aspace15_readdatavalid[0][0]),
      .avm_local_bb2_st__u2_inst0_writeack(local_avm_aspace15_writeack[0][0]),
      // AVM avm_local_bb4_ld__u10_inst0
      .avm_local_bb4_ld__u10_inst0_enable(local_avm_aspace15_enable[0][1]),
      .avm_local_bb4_ld__u10_inst0_read(local_avm_aspace15_read[0][1]),
      .avm_local_bb4_ld__u10_inst0_write(local_avm_aspace15_write[0][1]),
      .avm_local_bb4_ld__u10_inst0_burstcount(local_avm_aspace15_burstcount[0][1]),
      .avm_local_bb4_ld__u10_inst0_address(local_avm_aspace15_address[0][1]),
      .avm_local_bb4_ld__u10_inst0_writedata(local_avm_aspace15_writedata[0][1]),
      .avm_local_bb4_ld__u10_inst0_byteenable(local_avm_aspace15_byteenable[0][1]),
      .avm_local_bb4_ld__u10_inst0_waitrequest(local_avm_aspace15_waitrequest[0][1]),
      .avm_local_bb4_ld__u10_inst0_readdata(local_avm_aspace15_readdata[0][1]),
      .avm_local_bb4_ld__u10_inst0_readdatavalid(local_avm_aspace15_readdatavalid[0][1]),
      .avm_local_bb4_ld__u10_inst0_writeack(local_avm_aspace15_writeack[0][1]),
      // AVST avst_local_bb2__chan_Conf2Intere_active_inst0
      .avst_local_bb2__chan_Conf2Intere_active_inst0_valid(avm_channel_id_chan_Conf2Intere_active_read_valid),
      .avst_local_bb2__chan_Conf2Intere_active_inst0_ready(avm_channel_id_chan_Conf2Intere_active_read_ready),
      .avst_local_bb2__chan_Conf2Intere_active_inst0_data(avm_channel_id_chan_Conf2Intere_active_read_data),
      // AVST avst_local_bb2__chan_Conf2Intere_cnt_inst0
      .avst_local_bb2__chan_Conf2Intere_cnt_inst0_valid(avm_channel_id_chan_Conf2Intere_cnt_read_valid),
      .avst_local_bb2__chan_Conf2Intere_cnt_inst0_ready(avm_channel_id_chan_Conf2Intere_cnt_read_ready),
      .avst_local_bb2__chan_Conf2Intere_cnt_inst0_data(avm_channel_id_chan_Conf2Intere_cnt_read_data),
      // AVST avst_local_bb2__chan_Conf2Intere_mode_inst0
      .avst_local_bb2__chan_Conf2Intere_mode_inst0_valid(avm_channel_id_chan_Conf2Intere_mode_read_valid),
      .avst_local_bb2__chan_Conf2Intere_mode_inst0_ready(avm_channel_id_chan_Conf2Intere_mode_read_ready),
      .avst_local_bb2__chan_Conf2Intere_mode_inst0_data(avm_channel_id_chan_Conf2Intere_mode_read_data),
      // AVST avst_local_bb2__chan_Conf2Intere_x_inst0
      .avst_local_bb2__chan_Conf2Intere_x_inst0_valid(avm_channel_id_chan_Conf2Intere_x_read_valid),
      .avst_local_bb2__chan_Conf2Intere_x_inst0_ready(avm_channel_id_chan_Conf2Intere_x_read_ready),
      .avst_local_bb2__chan_Conf2Intere_x_inst0_data(avm_channel_id_chan_Conf2Intere_x_read_data),
      // AVST avst_local_bb2__chan_Conf2Intere_y_inst0
      .avst_local_bb2__chan_Conf2Intere_y_inst0_valid(avm_channel_id_chan_Conf2Intere_y_read_valid),
      .avst_local_bb2__chan_Conf2Intere_y_inst0_ready(avm_channel_id_chan_Conf2Intere_y_read_ready),
      .avst_local_bb2__chan_Conf2Intere_y_inst0_data(avm_channel_id_chan_Conf2Intere_y_read_data),
      // AVST avst_local_bb2__chan_Conf2Intere_z_inst0
      .avst_local_bb2__chan_Conf2Intere_z_inst0_valid(avm_channel_id_chan_Conf2Intere_z_read_valid),
      .avst_local_bb2__chan_Conf2Intere_z_inst0_ready(avm_channel_id_chan_Conf2Intere_z_read_ready),
      .avst_local_bb2__chan_Conf2Intere_z_inst0_data(avm_channel_id_chan_Conf2Intere_z_read_data),
      // AVST avst_local_bb6__chan_Intere2Store_active_inst0
      .avst_local_bb6__chan_Intere2Store_active_inst0_valid(avm_channel_id_chan_Intere2Store_active_write_valid),
      .avst_local_bb6__chan_Intere2Store_active_inst0_ready(avm_channel_id_chan_Intere2Store_active_write_ready),
      .avst_local_bb6__chan_Intere2Store_active_inst0_data(avm_channel_id_chan_Intere2Store_active_write_data),
      .avst_local_bb6__chan_Intere2Store_active_inst0_almostfull(avm_channel_id_chan_Intere2Store_active_write_almostfull),
      // AVST avst_local_bb6__chan_Intere2Store_cnt_inst0
      .avst_local_bb6__chan_Intere2Store_cnt_inst0_valid(avm_channel_id_chan_Intere2Store_cnt_write_valid),
      .avst_local_bb6__chan_Intere2Store_cnt_inst0_ready(avm_channel_id_chan_Intere2Store_cnt_write_ready),
      .avst_local_bb6__chan_Intere2Store_cnt_inst0_data(avm_channel_id_chan_Intere2Store_cnt_write_data),
      .avst_local_bb6__chan_Intere2Store_cnt_inst0_almostfull(avm_channel_id_chan_Intere2Store_cnt_write_almostfull),
      // AVST avst_local_bb6__chan_Intere2Store_intere_inst0
      .avst_local_bb6__chan_Intere2Store_intere_inst0_valid(avm_channel_id_chan_Intere2Store_intere_write_valid),
      .avst_local_bb6__chan_Intere2Store_intere_inst0_ready(avm_channel_id_chan_Intere2Store_intere_write_ready),
      .avst_local_bb6__chan_Intere2Store_intere_inst0_data(avm_channel_id_chan_Intere2Store_intere_write_data),
      .avst_local_bb6__chan_Intere2Store_intere_inst0_almostfull(avm_channel_id_chan_Intere2Store_intere_write_almostfull),
      // AVST avst_local_bb6__chan_Intere2Store_mode_inst0
      .avst_local_bb6__chan_Intere2Store_mode_inst0_valid(avm_channel_id_chan_Intere2Store_mode_write_valid),
      .avst_local_bb6__chan_Intere2Store_mode_inst0_ready(avm_channel_id_chan_Intere2Store_mode_write_ready),
      .avst_local_bb6__chan_Intere2Store_mode_inst0_data(avm_channel_id_chan_Intere2Store_mode_write_data),
      .avst_local_bb6__chan_Intere2Store_mode_inst0_almostfull(avm_channel_id_chan_Intere2Store_mode_write_almostfull),
      // AVM p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid)
   );

   assign lmem_invalid_single_bit = |lmem_invalid_aspaces;
   generate
   begin:local_mem_system_aspace11
      logic local_icm_arb_request [1][2];
      logic local_icm_arb_enable [1][2];
      logic local_icm_arb_read [1][2];
      logic local_icm_arb_write [1][2];
      logic local_icm_arb_burstcount [1][2];
      logic [6:0] local_icm_arb_address [1][2];
      logic [31:0] local_icm_arb_writedata [1][2];
      logic [3:0] local_icm_arb_byteenable [1][2];
      logic local_icm_arb_stall [1][2];
      logic local_icm_wrp_ack [1][2];
      logic local_icm_rrp_datavalid [1][2];
      logic [31:0] local_icm_rrp_data [1][2];

      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:local_mem_group
         for( __j = 0; __j < 2; __j = __j + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace11_enable[__i][__j]),
               .avm_read(local_avm_aspace11_read[__i][__j]),
               .avm_write(local_avm_aspace11_write[__i][__j]),
               .avm_burstcount(local_avm_aspace11_burstcount[__i][__j]),
               .avm_address(local_avm_aspace11_address[__i][__j]),
               .avm_writedata(local_avm_aspace11_writedata[__i][__j]),
               .avm_byteenable(local_avm_aspace11_byteenable[__i][__j]),
               .avm_waitrequest(local_avm_aspace11_waitrequest[__i][__j]),
               .avm_readdata(local_avm_aspace11_readdata[__i][__j]),
               .avm_readdatavalid(local_avm_aspace11_readdatavalid[__i][__j]),
               .avm_writeack(local_avm_aspace11_writeack[__i][__j]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__i][__j]),
               .ic_arb_enable(local_icm_arb_enable[__i][__j]),
               .ic_arb_read(local_icm_arb_read[__i][__j]),
               .ic_arb_write(local_icm_arb_write[__i][__j]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__i][__j]),
               .ic_arb_address(local_icm_arb_address[__i][__j]),
               .ic_arb_writedata(local_icm_arb_writedata[__i][__j]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__i][__j]),
               .ic_arb_stall(local_icm_arb_stall[__i][__j]),
               .ic_wrp_ack(local_icm_wrp_ack[__i][__j]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__i][__j]),
               .ic_rrp_data(local_icm_rrp_data[__i][__j])
            );

         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:bank
            logic port_enable [1:2];
            logic port_read [1:2];
            logic port_write [1:2];
            logic [6:0] port_address [1:2];
            logic [31:0] port_writedata [1:2];
            logic [3:0] port_byteenable [1:2];
            logic port_waitrequest [1:2];
            logic [31:0] port_readdata [1:2];
            logic port_readdatavalid [1:2];

            // INST mem0 of acl_mem1x
            acl_mem1x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .MEM_LATENCY(3),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("DUAL_PORT"),
               .PREFERRED_WIDTH(320),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2])
            );

         end

         for( __j = 0; __j < 2; __j = __j + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__i][__j]),
               .m_arb_enable(local_icm_arb_enable[__i][__j]),
               .m_arb_read(local_icm_arb_read[__i][__j]),
               .m_arb_write(local_icm_arb_write[__i][__j]),
               .m_arb_burstcount(local_icm_arb_burstcount[__i][__j]),
               .m_arb_address(local_icm_arb_address[__i][__j]),
               .m_arb_writedata(local_icm_arb_writedata[__i][__j]),
               .m_arb_byteenable(local_icm_arb_byteenable[__i][__j]),
               .m_arb_stall(local_icm_arb_stall[__i][__j]),
               .m_wrp_ack(local_icm_wrp_ack[__i][__j]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__i][__j]),
               .m_rrp_data(local_icm_rrp_data[__i][__j]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_4
            Krnl_GA_system_interconnect_4 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_5
            Krnl_GA_system_interconnect_5 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

      end

   end
   endgenerate

   generate
   begin:local_mem_system_aspace13
      logic local_icm_arb_request [1][2];
      logic local_icm_arb_enable [1][2];
      logic local_icm_arb_read [1][2];
      logic local_icm_arb_write [1][2];
      logic local_icm_arb_burstcount [1][2];
      logic [6:0] local_icm_arb_address [1][2];
      logic [31:0] local_icm_arb_writedata [1][2];
      logic [3:0] local_icm_arb_byteenable [1][2];
      logic local_icm_arb_stall [1][2];
      logic local_icm_wrp_ack [1][2];
      logic local_icm_rrp_datavalid [1][2];
      logic [31:0] local_icm_rrp_data [1][2];

      for( __j = 0; __j < 1; __j = __j + 1 )
      begin:local_mem_group
         for( __k = 0; __k < 2; __k = __k + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace13_enable[__j][__k]),
               .avm_read(local_avm_aspace13_read[__j][__k]),
               .avm_write(local_avm_aspace13_write[__j][__k]),
               .avm_burstcount(local_avm_aspace13_burstcount[__j][__k]),
               .avm_address(local_avm_aspace13_address[__j][__k]),
               .avm_writedata(local_avm_aspace13_writedata[__j][__k]),
               .avm_byteenable(local_avm_aspace13_byteenable[__j][__k]),
               .avm_waitrequest(local_avm_aspace13_waitrequest[__j][__k]),
               .avm_readdata(local_avm_aspace13_readdata[__j][__k]),
               .avm_readdatavalid(local_avm_aspace13_readdatavalid[__j][__k]),
               .avm_writeack(local_avm_aspace13_writeack[__j][__k]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__j][__k]),
               .ic_arb_enable(local_icm_arb_enable[__j][__k]),
               .ic_arb_read(local_icm_arb_read[__j][__k]),
               .ic_arb_write(local_icm_arb_write[__j][__k]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__j][__k]),
               .ic_arb_address(local_icm_arb_address[__j][__k]),
               .ic_arb_writedata(local_icm_arb_writedata[__j][__k]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__j][__k]),
               .ic_arb_stall(local_icm_arb_stall[__j][__k]),
               .ic_wrp_ack(local_icm_wrp_ack[__j][__k]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__j][__k]),
               .ic_rrp_data(local_icm_rrp_data[__j][__k])
            );

         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:bank
            logic port_enable [1:2];
            logic port_read [1:2];
            logic port_write [1:2];
            logic [6:0] port_address [1:2];
            logic [31:0] port_writedata [1:2];
            logic [3:0] port_byteenable [1:2];
            logic port_waitrequest [1:2];
            logic [31:0] port_readdata [1:2];
            logic port_readdatavalid [1:2];

            // INST mem0 of acl_mem1x
            acl_mem1x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .MEM_LATENCY(3),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("DUAL_PORT"),
               .PREFERRED_WIDTH(320),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2])
            );

         end

         for( __k = 0; __k < 2; __k = __k + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__j][__k]),
               .m_arb_enable(local_icm_arb_enable[__j][__k]),
               .m_arb_read(local_icm_arb_read[__j][__k]),
               .m_arb_write(local_icm_arb_write[__j][__k]),
               .m_arb_burstcount(local_icm_arb_burstcount[__j][__k]),
               .m_arb_address(local_icm_arb_address[__j][__k]),
               .m_arb_writedata(local_icm_arb_writedata[__j][__k]),
               .m_arb_byteenable(local_icm_arb_byteenable[__j][__k]),
               .m_arb_stall(local_icm_arb_stall[__j][__k]),
               .m_wrp_ack(local_icm_wrp_ack[__j][__k]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__j][__k]),
               .m_rrp_data(local_icm_rrp_data[__j][__k]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_4
            Krnl_GA_system_interconnect_4 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_5
            Krnl_GA_system_interconnect_5 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

      end

   end
   endgenerate

   generate
   begin:local_mem_system_aspace15
      logic local_icm_arb_request [1][2];
      logic local_icm_arb_enable [1][2];
      logic local_icm_arb_read [1][2];
      logic local_icm_arb_write [1][2];
      logic local_icm_arb_burstcount [1][2];
      logic [6:0] local_icm_arb_address [1][2];
      logic [31:0] local_icm_arb_writedata [1][2];
      logic [3:0] local_icm_arb_byteenable [1][2];
      logic local_icm_arb_stall [1][2];
      logic local_icm_wrp_ack [1][2];
      logic local_icm_rrp_datavalid [1][2];
      logic [31:0] local_icm_rrp_data [1][2];

      for( __k = 0; __k < 1; __k = __k + 1 )
      begin:local_mem_group
         for( __l = 0; __l < 2; __l = __l + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace15_enable[__k][__l]),
               .avm_read(local_avm_aspace15_read[__k][__l]),
               .avm_write(local_avm_aspace15_write[__k][__l]),
               .avm_burstcount(local_avm_aspace15_burstcount[__k][__l]),
               .avm_address(local_avm_aspace15_address[__k][__l]),
               .avm_writedata(local_avm_aspace15_writedata[__k][__l]),
               .avm_byteenable(local_avm_aspace15_byteenable[__k][__l]),
               .avm_waitrequest(local_avm_aspace15_waitrequest[__k][__l]),
               .avm_readdata(local_avm_aspace15_readdata[__k][__l]),
               .avm_readdatavalid(local_avm_aspace15_readdatavalid[__k][__l]),
               .avm_writeack(local_avm_aspace15_writeack[__k][__l]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__k][__l]),
               .ic_arb_enable(local_icm_arb_enable[__k][__l]),
               .ic_arb_read(local_icm_arb_read[__k][__l]),
               .ic_arb_write(local_icm_arb_write[__k][__l]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__k][__l]),
               .ic_arb_address(local_icm_arb_address[__k][__l]),
               .ic_arb_writedata(local_icm_arb_writedata[__k][__l]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__k][__l]),
               .ic_arb_stall(local_icm_arb_stall[__k][__l]),
               .ic_wrp_ack(local_icm_wrp_ack[__k][__l]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__k][__l]),
               .ic_rrp_data(local_icm_rrp_data[__k][__l])
            );

         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:bank
            logic port_enable [1:2];
            logic port_read [1:2];
            logic port_write [1:2];
            logic [6:0] port_address [1:2];
            logic [31:0] port_writedata [1:2];
            logic [3:0] port_byteenable [1:2];
            logic port_waitrequest [1:2];
            logic [31:0] port_readdata [1:2];
            logic port_readdatavalid [1:2];

            // INST mem0 of acl_mem1x
            acl_mem1x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .MEM_LATENCY(3),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("DUAL_PORT"),
               .PREFERRED_WIDTH(320),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2])
            );

         end

         for( __l = 0; __l < 2; __l = __l + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__k][__l]),
               .m_arb_enable(local_icm_arb_enable[__k][__l]),
               .m_arb_read(local_icm_arb_read[__k][__l]),
               .m_arb_write(local_icm_arb_write[__k][__l]),
               .m_arb_burstcount(local_icm_arb_burstcount[__k][__l]),
               .m_arb_address(local_icm_arb_address[__k][__l]),
               .m_arb_writedata(local_icm_arb_writedata[__k][__l]),
               .m_arb_byteenable(local_icm_arb_byteenable[__k][__l]),
               .m_arb_stall(local_icm_arb_stall[__k][__l]),
               .m_wrp_ack(local_icm_wrp_ack[__k][__l]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__k][__l]),
               .m_rrp_data(local_icm_rrp_data[__k][__l]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_4
            Krnl_GA_system_interconnect_4 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_5
            Krnl_GA_system_interconnect_5 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

      end

   end
   endgenerate

endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_IntraE_top_wrapper_0
/////////////////////////////////////////////////////////////////
module Krnl_IntraE_top_wrapper_0
(
   input logic start,
   input logic [223:0] kernel_arguments,
   input logic [31:0] work_dim,
   input logic [31:0] global_offset [2:0],
   output logic kernel_valid_out,
   output logic has_a_write_pending,
   output logic has_a_lsu_active,
   input logic [31:0] global_id [2:0],
   input logic [31:0] local_id [2:0],
   input logic [31:0] group_id [2:0],
   input logic [31:0] global_size [2:0],
   input logic [31:0] local_size [2:0],
   input logic [31:0] num_groups [2:0],
   input logic [31:0] workgroup_size,
   output logic kernel_stall_out,
   input logic kernel_valid_in,
   input logic clock,
   input logic resetn,
   input logic clock2x,
   // AVM avm_local_bb0_ld__inst0
   output logic avm_local_bb0_ld__inst0_enable,
   output logic avm_local_bb0_ld__inst0_read,
   output logic avm_local_bb0_ld__inst0_write,
   output logic [4:0] avm_local_bb0_ld__inst0_burstcount,
   output logic [30:0] avm_local_bb0_ld__inst0_address,
   output logic [511:0] avm_local_bb0_ld__inst0_writedata,
   output logic [63:0] avm_local_bb0_ld__inst0_byteenable,
   input logic avm_local_bb0_ld__inst0_waitrequest,
   input logic [511:0] avm_local_bb0_ld__inst0_readdata,
   input logic avm_local_bb0_ld__inst0_readdatavalid,
   input logic avm_local_bb0_ld__inst0_writeack,
   // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable,
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read,
   output logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write,
   output logic [4:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount,
   output logic [30:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address,
   output logic [511:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata,
   output logic [63:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest,
   input logic [511:0] avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid,
   input logic avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack,
   // AVM avm_local_bb3_st__inst0
   output logic avm_local_bb3_st__inst0_enable,
   output logic avm_local_bb3_st__inst0_read,
   output logic avm_local_bb3_st__inst0_write,
   output logic [4:0] avm_local_bb3_st__inst0_burstcount,
   output logic [30:0] avm_local_bb3_st__inst0_address,
   output logic [511:0] avm_local_bb3_st__inst0_writedata,
   output logic [63:0] avm_local_bb3_st__inst0_byteenable,
   input logic avm_local_bb3_st__inst0_waitrequest,
   input logic [511:0] avm_local_bb3_st__inst0_readdata,
   input logic avm_local_bb3_st__inst0_readdatavalid,
   input logic avm_local_bb3_st__inst0_writeack,
   // AVM avm_local_bb4_ld__u18_inst0
   output logic avm_local_bb4_ld__u18_inst0_enable,
   output logic avm_local_bb4_ld__u18_inst0_read,
   output logic avm_local_bb4_ld__u18_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u18_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u18_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u18_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u18_inst0_byteenable,
   input logic avm_local_bb4_ld__u18_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u18_inst0_readdata,
   input logic avm_local_bb4_ld__u18_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u18_inst0_writeack,
   // AVM avm_local_bb4_ld__u19_inst0
   output logic avm_local_bb4_ld__u19_inst0_enable,
   output logic avm_local_bb4_ld__u19_inst0_read,
   output logic avm_local_bb4_ld__u19_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u19_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u19_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u19_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u19_inst0_byteenable,
   input logic avm_local_bb4_ld__u19_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u19_inst0_readdata,
   input logic avm_local_bb4_ld__u19_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u19_inst0_writeack,
   // AVM avm_local_bb4_ld__u20_inst0
   output logic avm_local_bb4_ld__u20_inst0_enable,
   output logic avm_local_bb4_ld__u20_inst0_read,
   output logic avm_local_bb4_ld__u20_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u20_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u20_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u20_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u20_inst0_byteenable,
   input logic avm_local_bb4_ld__u20_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u20_inst0_readdata,
   input logic avm_local_bb4_ld__u20_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u20_inst0_writeack,
   // AVM avm_local_bb4_ld__u21_inst0
   output logic avm_local_bb4_ld__u21_inst0_enable,
   output logic avm_local_bb4_ld__u21_inst0_read,
   output logic avm_local_bb4_ld__u21_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u21_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u21_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u21_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u21_inst0_byteenable,
   input logic avm_local_bb4_ld__u21_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u21_inst0_readdata,
   input logic avm_local_bb4_ld__u21_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u21_inst0_writeack,
   // AVM avm_local_bb4_ld__u22_inst0
   output logic avm_local_bb4_ld__u22_inst0_enable,
   output logic avm_local_bb4_ld__u22_inst0_read,
   output logic avm_local_bb4_ld__u22_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u22_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u22_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u22_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u22_inst0_byteenable,
   input logic avm_local_bb4_ld__u22_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u22_inst0_readdata,
   input logic avm_local_bb4_ld__u22_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u22_inst0_writeack,
   // AVM avm_local_bb4_ld__u23_inst0
   output logic avm_local_bb4_ld__u23_inst0_enable,
   output logic avm_local_bb4_ld__u23_inst0_read,
   output logic avm_local_bb4_ld__u23_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u23_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u23_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u23_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u23_inst0_byteenable,
   input logic avm_local_bb4_ld__u23_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u23_inst0_readdata,
   input logic avm_local_bb4_ld__u23_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u23_inst0_writeack,
   // AVM avm_local_bb4_ld__u24_inst0
   output logic avm_local_bb4_ld__u24_inst0_enable,
   output logic avm_local_bb4_ld__u24_inst0_read,
   output logic avm_local_bb4_ld__u24_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u24_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u24_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u24_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u24_inst0_byteenable,
   input logic avm_local_bb4_ld__u24_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u24_inst0_readdata,
   input logic avm_local_bb4_ld__u24_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u24_inst0_writeack,
   // AVM avm_local_bb4_ld__u25_inst0
   output logic avm_local_bb4_ld__u25_inst0_enable,
   output logic avm_local_bb4_ld__u25_inst0_read,
   output logic avm_local_bb4_ld__u25_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u25_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u25_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u25_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u25_inst0_byteenable,
   input logic avm_local_bb4_ld__u25_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u25_inst0_readdata,
   input logic avm_local_bb4_ld__u25_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u25_inst0_writeack,
   // AVM avm_local_bb4_ld__u26_inst0
   output logic avm_local_bb4_ld__u26_inst0_enable,
   output logic avm_local_bb4_ld__u26_inst0_read,
   output logic avm_local_bb4_ld__u26_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u26_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u26_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u26_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u26_inst0_byteenable,
   input logic avm_local_bb4_ld__u26_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u26_inst0_readdata,
   input logic avm_local_bb4_ld__u26_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u26_inst0_writeack,
   // AVM avm_local_bb4_ld__u27_inst0
   output logic avm_local_bb4_ld__u27_inst0_enable,
   output logic avm_local_bb4_ld__u27_inst0_read,
   output logic avm_local_bb4_ld__u27_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u27_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u27_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u27_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u27_inst0_byteenable,
   input logic avm_local_bb4_ld__u27_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u27_inst0_readdata,
   input logic avm_local_bb4_ld__u27_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u27_inst0_writeack,
   // AVM avm_local_bb4_ld__u28_inst0
   output logic avm_local_bb4_ld__u28_inst0_enable,
   output logic avm_local_bb4_ld__u28_inst0_read,
   output logic avm_local_bb4_ld__u28_inst0_write,
   output logic [4:0] avm_local_bb4_ld__u28_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld__u28_inst0_address,
   output logic [511:0] avm_local_bb4_ld__u28_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld__u28_inst0_byteenable,
   input logic avm_local_bb4_ld__u28_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld__u28_inst0_readdata,
   input logic avm_local_bb4_ld__u28_inst0_readdatavalid,
   input logic avm_local_bb4_ld__u28_inst0_writeack,
   // AVM avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0
   output logic avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_enable,
   output logic avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_read,
   output logic avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_write,
   output logic [4:0] avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_burstcount,
   output logic [30:0] avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_address,
   output logic [511:0] avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writedata,
   output logic [63:0] avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_byteenable,
   input logic avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_waitrequest,
   input logic [511:0] avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdata,
   input logic avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid,
   input logic avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writeack,
   // AVST avm_channel_id_chan_Conf2Intrae_active_read
   input logic avm_channel_id_chan_Conf2Intrae_active_read_valid,
   output logic avm_channel_id_chan_Conf2Intrae_active_read_ready,
   input logic [7:0] avm_channel_id_chan_Conf2Intrae_active_read_data,
   // AVST avm_channel_id_chan_Conf2Intrae_cnt_read
   input logic avm_channel_id_chan_Conf2Intrae_cnt_read_valid,
   output logic avm_channel_id_chan_Conf2Intrae_cnt_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intrae_cnt_read_data,
   // AVST avm_channel_id_chan_Conf2Intrae_mode_read
   input logic avm_channel_id_chan_Conf2Intrae_mode_read_valid,
   output logic avm_channel_id_chan_Conf2Intrae_mode_read_ready,
   input logic [7:0] avm_channel_id_chan_Conf2Intrae_mode_read_data,
   // AVST avm_channel_id_chan_Conf2Intrae_x_read
   input logic avm_channel_id_chan_Conf2Intrae_x_read_valid,
   output logic avm_channel_id_chan_Conf2Intrae_x_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intrae_x_read_data,
   // AVST avm_channel_id_chan_Conf2Intrae_y_read
   input logic avm_channel_id_chan_Conf2Intrae_y_read_valid,
   output logic avm_channel_id_chan_Conf2Intrae_y_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intrae_y_read_data,
   // AVST avm_channel_id_chan_Conf2Intrae_z_read
   input logic avm_channel_id_chan_Conf2Intrae_z_read_valid,
   output logic avm_channel_id_chan_Conf2Intrae_z_read_ready,
   input logic [31:0] avm_channel_id_chan_Conf2Intrae_z_read_data,
   // AVST avm_channel_id_chan_Intrae2Store_active_write
   output logic avm_channel_id_chan_Intrae2Store_active_write_valid,
   input logic avm_channel_id_chan_Intrae2Store_active_write_ready,
   output logic [7:0] avm_channel_id_chan_Intrae2Store_active_write_data,
   input logic avm_channel_id_chan_Intrae2Store_active_write_almostfull,
   // AVST avm_channel_id_chan_Intrae2Store_cnt_write
   output logic avm_channel_id_chan_Intrae2Store_cnt_write_valid,
   input logic avm_channel_id_chan_Intrae2Store_cnt_write_ready,
   output logic [31:0] avm_channel_id_chan_Intrae2Store_cnt_write_data,
   input logic avm_channel_id_chan_Intrae2Store_cnt_write_almostfull,
   // AVST avm_channel_id_chan_Intrae2Store_intrae_write
   output logic avm_channel_id_chan_Intrae2Store_intrae_write_valid,
   input logic avm_channel_id_chan_Intrae2Store_intrae_write_ready,
   output logic [31:0] avm_channel_id_chan_Intrae2Store_intrae_write_data,
   input logic avm_channel_id_chan_Intrae2Store_intrae_write_almostfull,
   // AVST avm_channel_id_chan_Intrae2Store_mode_write
   output logic avm_channel_id_chan_Intrae2Store_mode_write_valid,
   input logic avm_channel_id_chan_Intrae2Store_mode_write_ready,
   output logic [7:0] avm_channel_id_chan_Intrae2Store_mode_write_data,
   input logic avm_channel_id_chan_Intrae2Store_mode_write_almostfull,
   // AVM p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid
);
   genvar __i;
   genvar __j;
   genvar __k;
   genvar __l;
   logic lmem_invalid_single_bit;
   logic [2:0] lmem_invalid_aspaces;
   logic local_avm_aspace17_enable [1][3];
   logic local_avm_aspace17_read [1][3];
   logic local_avm_aspace17_write [1][3];
   logic local_avm_aspace17_burstcount [1][3];
   logic [31:0] local_avm_aspace17_address [1][3];
   logic [31:0] local_avm_aspace17_writedata [1][3];
   logic [3:0] local_avm_aspace17_byteenable [1][3];
   logic local_avm_aspace17_waitrequest [1][3];
   logic [31:0] local_avm_aspace17_readdata [1][3];
   logic local_avm_aspace17_readdatavalid [1][3];
   logic local_avm_aspace17_writeack [1][3];
   logic local_avm_aspace19_enable [1][3];
   logic local_avm_aspace19_read [1][3];
   logic local_avm_aspace19_write [1][3];
   logic local_avm_aspace19_burstcount [1][3];
   logic [31:0] local_avm_aspace19_address [1][3];
   logic [31:0] local_avm_aspace19_writedata [1][3];
   logic [3:0] local_avm_aspace19_byteenable [1][3];
   logic local_avm_aspace19_waitrequest [1][3];
   logic [31:0] local_avm_aspace19_readdata [1][3];
   logic local_avm_aspace19_readdatavalid [1][3];
   logic local_avm_aspace19_writeack [1][3];
   logic local_avm_aspace21_enable [1][3];
   logic local_avm_aspace21_read [1][3];
   logic local_avm_aspace21_write [1][3];
   logic local_avm_aspace21_burstcount [1][3];
   logic [31:0] local_avm_aspace21_address [1][3];
   logic [31:0] local_avm_aspace21_writedata [1][3];
   logic [3:0] local_avm_aspace21_byteenable [1][3];
   logic local_avm_aspace21_waitrequest [1][3];
   logic [31:0] local_avm_aspace21_readdata [1][3];
   logic local_avm_aspace21_readdatavalid [1][3];
   logic local_avm_aspace21_writeack [1][3];

   // INST kernel of Krnl_IntraE_function_wrapper
   Krnl_IntraE_function_wrapper kernel
   (
      .local_router_hang(lmem_invalid_single_bit),
      .start(start),
      .kernel_arguments(kernel_arguments),
      .work_dim(work_dim),
      .global_offset_0(global_offset[0]),
      .global_offset_1(global_offset[1]),
      .global_offset_2(global_offset[2]),
      .kernel_valid_out(kernel_valid_out),
      .has_a_write_pending(has_a_write_pending),
      .has_a_lsu_active(has_a_lsu_active),
      .global_id_0(global_id[0]),
      .global_id_1(global_id[1]),
      .global_id_2(global_id[2]),
      .local_id_0(local_id[0]),
      .local_id_1(local_id[1]),
      .local_id_2(local_id[2]),
      .group_id_0(group_id[0]),
      .group_id_1(group_id[1]),
      .group_id_2(group_id[2]),
      .global_size_0(global_size[0]),
      .global_size_1(global_size[1]),
      .global_size_2(global_size[2]),
      .local_size_0(local_size[0]),
      .local_size_1(local_size[1]),
      .local_size_2(local_size[2]),
      .num_groups_0(num_groups[0]),
      .num_groups_1(num_groups[1]),
      .num_groups_2(num_groups[2]),
      .workgroup_size(workgroup_size),
      .kernel_stall_out(kernel_stall_out),
      .kernel_valid_in(kernel_valid_in),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb0_ld__inst0
      .avm_local_bb0_ld__inst0_enable(avm_local_bb0_ld__inst0_enable),
      .avm_local_bb0_ld__inst0_read(avm_local_bb0_ld__inst0_read),
      .avm_local_bb0_ld__inst0_write(avm_local_bb0_ld__inst0_write),
      .avm_local_bb0_ld__inst0_burstcount(avm_local_bb0_ld__inst0_burstcount),
      .avm_local_bb0_ld__inst0_address(avm_local_bb0_ld__inst0_address),
      .avm_local_bb0_ld__inst0_writedata(avm_local_bb0_ld__inst0_writedata),
      .avm_local_bb0_ld__inst0_byteenable(avm_local_bb0_ld__inst0_byteenable),
      .avm_local_bb0_ld__inst0_waitrequest(avm_local_bb0_ld__inst0_waitrequest),
      .avm_local_bb0_ld__inst0_readdata(avm_local_bb0_ld__inst0_readdata),
      .avm_local_bb0_ld__inst0_readdatavalid(avm_local_bb0_ld__inst0_readdatavalid),
      .avm_local_bb0_ld__inst0_writeack(avm_local_bb0_ld__inst0_writeack),
      // AVM avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_enable),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_read),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_write),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_burstcount),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_address),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writedata),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_byteenable),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_waitrequest),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdata),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_readdatavalid),
      .avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack(avm_local_bb0_ld_memcoalesce_DockConst_load_0_inst0_writeack),
      // AVM avm_local_bb3_st__inst0
      .avm_local_bb3_st__inst0_enable(avm_local_bb3_st__inst0_enable),
      .avm_local_bb3_st__inst0_read(avm_local_bb3_st__inst0_read),
      .avm_local_bb3_st__inst0_write(avm_local_bb3_st__inst0_write),
      .avm_local_bb3_st__inst0_burstcount(avm_local_bb3_st__inst0_burstcount),
      .avm_local_bb3_st__inst0_address(avm_local_bb3_st__inst0_address),
      .avm_local_bb3_st__inst0_writedata(avm_local_bb3_st__inst0_writedata),
      .avm_local_bb3_st__inst0_byteenable(avm_local_bb3_st__inst0_byteenable),
      .avm_local_bb3_st__inst0_waitrequest(avm_local_bb3_st__inst0_waitrequest),
      .avm_local_bb3_st__inst0_readdata(avm_local_bb3_st__inst0_readdata),
      .avm_local_bb3_st__inst0_readdatavalid(avm_local_bb3_st__inst0_readdatavalid),
      .avm_local_bb3_st__inst0_writeack(avm_local_bb3_st__inst0_writeack),
      // AVM avm_local_bb4_ld__u18_inst0
      .avm_local_bb4_ld__u18_inst0_enable(avm_local_bb4_ld__u18_inst0_enable),
      .avm_local_bb4_ld__u18_inst0_read(avm_local_bb4_ld__u18_inst0_read),
      .avm_local_bb4_ld__u18_inst0_write(avm_local_bb4_ld__u18_inst0_write),
      .avm_local_bb4_ld__u18_inst0_burstcount(avm_local_bb4_ld__u18_inst0_burstcount),
      .avm_local_bb4_ld__u18_inst0_address(avm_local_bb4_ld__u18_inst0_address),
      .avm_local_bb4_ld__u18_inst0_writedata(avm_local_bb4_ld__u18_inst0_writedata),
      .avm_local_bb4_ld__u18_inst0_byteenable(avm_local_bb4_ld__u18_inst0_byteenable),
      .avm_local_bb4_ld__u18_inst0_waitrequest(avm_local_bb4_ld__u18_inst0_waitrequest),
      .avm_local_bb4_ld__u18_inst0_readdata(avm_local_bb4_ld__u18_inst0_readdata),
      .avm_local_bb4_ld__u18_inst0_readdatavalid(avm_local_bb4_ld__u18_inst0_readdatavalid),
      .avm_local_bb4_ld__u18_inst0_writeack(avm_local_bb4_ld__u18_inst0_writeack),
      // AVM avm_local_bb4_ld__u19_inst0
      .avm_local_bb4_ld__u19_inst0_enable(avm_local_bb4_ld__u19_inst0_enable),
      .avm_local_bb4_ld__u19_inst0_read(avm_local_bb4_ld__u19_inst0_read),
      .avm_local_bb4_ld__u19_inst0_write(avm_local_bb4_ld__u19_inst0_write),
      .avm_local_bb4_ld__u19_inst0_burstcount(avm_local_bb4_ld__u19_inst0_burstcount),
      .avm_local_bb4_ld__u19_inst0_address(avm_local_bb4_ld__u19_inst0_address),
      .avm_local_bb4_ld__u19_inst0_writedata(avm_local_bb4_ld__u19_inst0_writedata),
      .avm_local_bb4_ld__u19_inst0_byteenable(avm_local_bb4_ld__u19_inst0_byteenable),
      .avm_local_bb4_ld__u19_inst0_waitrequest(avm_local_bb4_ld__u19_inst0_waitrequest),
      .avm_local_bb4_ld__u19_inst0_readdata(avm_local_bb4_ld__u19_inst0_readdata),
      .avm_local_bb4_ld__u19_inst0_readdatavalid(avm_local_bb4_ld__u19_inst0_readdatavalid),
      .avm_local_bb4_ld__u19_inst0_writeack(avm_local_bb4_ld__u19_inst0_writeack),
      // AVM avm_local_bb4_ld__u20_inst0
      .avm_local_bb4_ld__u20_inst0_enable(avm_local_bb4_ld__u20_inst0_enable),
      .avm_local_bb4_ld__u20_inst0_read(avm_local_bb4_ld__u20_inst0_read),
      .avm_local_bb4_ld__u20_inst0_write(avm_local_bb4_ld__u20_inst0_write),
      .avm_local_bb4_ld__u20_inst0_burstcount(avm_local_bb4_ld__u20_inst0_burstcount),
      .avm_local_bb4_ld__u20_inst0_address(avm_local_bb4_ld__u20_inst0_address),
      .avm_local_bb4_ld__u20_inst0_writedata(avm_local_bb4_ld__u20_inst0_writedata),
      .avm_local_bb4_ld__u20_inst0_byteenable(avm_local_bb4_ld__u20_inst0_byteenable),
      .avm_local_bb4_ld__u20_inst0_waitrequest(avm_local_bb4_ld__u20_inst0_waitrequest),
      .avm_local_bb4_ld__u20_inst0_readdata(avm_local_bb4_ld__u20_inst0_readdata),
      .avm_local_bb4_ld__u20_inst0_readdatavalid(avm_local_bb4_ld__u20_inst0_readdatavalid),
      .avm_local_bb4_ld__u20_inst0_writeack(avm_local_bb4_ld__u20_inst0_writeack),
      // AVM avm_local_bb4_ld__u21_inst0
      .avm_local_bb4_ld__u21_inst0_enable(avm_local_bb4_ld__u21_inst0_enable),
      .avm_local_bb4_ld__u21_inst0_read(avm_local_bb4_ld__u21_inst0_read),
      .avm_local_bb4_ld__u21_inst0_write(avm_local_bb4_ld__u21_inst0_write),
      .avm_local_bb4_ld__u21_inst0_burstcount(avm_local_bb4_ld__u21_inst0_burstcount),
      .avm_local_bb4_ld__u21_inst0_address(avm_local_bb4_ld__u21_inst0_address),
      .avm_local_bb4_ld__u21_inst0_writedata(avm_local_bb4_ld__u21_inst0_writedata),
      .avm_local_bb4_ld__u21_inst0_byteenable(avm_local_bb4_ld__u21_inst0_byteenable),
      .avm_local_bb4_ld__u21_inst0_waitrequest(avm_local_bb4_ld__u21_inst0_waitrequest),
      .avm_local_bb4_ld__u21_inst0_readdata(avm_local_bb4_ld__u21_inst0_readdata),
      .avm_local_bb4_ld__u21_inst0_readdatavalid(avm_local_bb4_ld__u21_inst0_readdatavalid),
      .avm_local_bb4_ld__u21_inst0_writeack(avm_local_bb4_ld__u21_inst0_writeack),
      // AVM avm_local_bb4_ld__u22_inst0
      .avm_local_bb4_ld__u22_inst0_enable(avm_local_bb4_ld__u22_inst0_enable),
      .avm_local_bb4_ld__u22_inst0_read(avm_local_bb4_ld__u22_inst0_read),
      .avm_local_bb4_ld__u22_inst0_write(avm_local_bb4_ld__u22_inst0_write),
      .avm_local_bb4_ld__u22_inst0_burstcount(avm_local_bb4_ld__u22_inst0_burstcount),
      .avm_local_bb4_ld__u22_inst0_address(avm_local_bb4_ld__u22_inst0_address),
      .avm_local_bb4_ld__u22_inst0_writedata(avm_local_bb4_ld__u22_inst0_writedata),
      .avm_local_bb4_ld__u22_inst0_byteenable(avm_local_bb4_ld__u22_inst0_byteenable),
      .avm_local_bb4_ld__u22_inst0_waitrequest(avm_local_bb4_ld__u22_inst0_waitrequest),
      .avm_local_bb4_ld__u22_inst0_readdata(avm_local_bb4_ld__u22_inst0_readdata),
      .avm_local_bb4_ld__u22_inst0_readdatavalid(avm_local_bb4_ld__u22_inst0_readdatavalid),
      .avm_local_bb4_ld__u22_inst0_writeack(avm_local_bb4_ld__u22_inst0_writeack),
      // AVM avm_local_bb4_ld__u23_inst0
      .avm_local_bb4_ld__u23_inst0_enable(avm_local_bb4_ld__u23_inst0_enable),
      .avm_local_bb4_ld__u23_inst0_read(avm_local_bb4_ld__u23_inst0_read),
      .avm_local_bb4_ld__u23_inst0_write(avm_local_bb4_ld__u23_inst0_write),
      .avm_local_bb4_ld__u23_inst0_burstcount(avm_local_bb4_ld__u23_inst0_burstcount),
      .avm_local_bb4_ld__u23_inst0_address(avm_local_bb4_ld__u23_inst0_address),
      .avm_local_bb4_ld__u23_inst0_writedata(avm_local_bb4_ld__u23_inst0_writedata),
      .avm_local_bb4_ld__u23_inst0_byteenable(avm_local_bb4_ld__u23_inst0_byteenable),
      .avm_local_bb4_ld__u23_inst0_waitrequest(avm_local_bb4_ld__u23_inst0_waitrequest),
      .avm_local_bb4_ld__u23_inst0_readdata(avm_local_bb4_ld__u23_inst0_readdata),
      .avm_local_bb4_ld__u23_inst0_readdatavalid(avm_local_bb4_ld__u23_inst0_readdatavalid),
      .avm_local_bb4_ld__u23_inst0_writeack(avm_local_bb4_ld__u23_inst0_writeack),
      // AVM avm_local_bb4_ld__u24_inst0
      .avm_local_bb4_ld__u24_inst0_enable(avm_local_bb4_ld__u24_inst0_enable),
      .avm_local_bb4_ld__u24_inst0_read(avm_local_bb4_ld__u24_inst0_read),
      .avm_local_bb4_ld__u24_inst0_write(avm_local_bb4_ld__u24_inst0_write),
      .avm_local_bb4_ld__u24_inst0_burstcount(avm_local_bb4_ld__u24_inst0_burstcount),
      .avm_local_bb4_ld__u24_inst0_address(avm_local_bb4_ld__u24_inst0_address),
      .avm_local_bb4_ld__u24_inst0_writedata(avm_local_bb4_ld__u24_inst0_writedata),
      .avm_local_bb4_ld__u24_inst0_byteenable(avm_local_bb4_ld__u24_inst0_byteenable),
      .avm_local_bb4_ld__u24_inst0_waitrequest(avm_local_bb4_ld__u24_inst0_waitrequest),
      .avm_local_bb4_ld__u24_inst0_readdata(avm_local_bb4_ld__u24_inst0_readdata),
      .avm_local_bb4_ld__u24_inst0_readdatavalid(avm_local_bb4_ld__u24_inst0_readdatavalid),
      .avm_local_bb4_ld__u24_inst0_writeack(avm_local_bb4_ld__u24_inst0_writeack),
      // AVM avm_local_bb4_ld__u25_inst0
      .avm_local_bb4_ld__u25_inst0_enable(avm_local_bb4_ld__u25_inst0_enable),
      .avm_local_bb4_ld__u25_inst0_read(avm_local_bb4_ld__u25_inst0_read),
      .avm_local_bb4_ld__u25_inst0_write(avm_local_bb4_ld__u25_inst0_write),
      .avm_local_bb4_ld__u25_inst0_burstcount(avm_local_bb4_ld__u25_inst0_burstcount),
      .avm_local_bb4_ld__u25_inst0_address(avm_local_bb4_ld__u25_inst0_address),
      .avm_local_bb4_ld__u25_inst0_writedata(avm_local_bb4_ld__u25_inst0_writedata),
      .avm_local_bb4_ld__u25_inst0_byteenable(avm_local_bb4_ld__u25_inst0_byteenable),
      .avm_local_bb4_ld__u25_inst0_waitrequest(avm_local_bb4_ld__u25_inst0_waitrequest),
      .avm_local_bb4_ld__u25_inst0_readdata(avm_local_bb4_ld__u25_inst0_readdata),
      .avm_local_bb4_ld__u25_inst0_readdatavalid(avm_local_bb4_ld__u25_inst0_readdatavalid),
      .avm_local_bb4_ld__u25_inst0_writeack(avm_local_bb4_ld__u25_inst0_writeack),
      // AVM avm_local_bb4_ld__u26_inst0
      .avm_local_bb4_ld__u26_inst0_enable(avm_local_bb4_ld__u26_inst0_enable),
      .avm_local_bb4_ld__u26_inst0_read(avm_local_bb4_ld__u26_inst0_read),
      .avm_local_bb4_ld__u26_inst0_write(avm_local_bb4_ld__u26_inst0_write),
      .avm_local_bb4_ld__u26_inst0_burstcount(avm_local_bb4_ld__u26_inst0_burstcount),
      .avm_local_bb4_ld__u26_inst0_address(avm_local_bb4_ld__u26_inst0_address),
      .avm_local_bb4_ld__u26_inst0_writedata(avm_local_bb4_ld__u26_inst0_writedata),
      .avm_local_bb4_ld__u26_inst0_byteenable(avm_local_bb4_ld__u26_inst0_byteenable),
      .avm_local_bb4_ld__u26_inst0_waitrequest(avm_local_bb4_ld__u26_inst0_waitrequest),
      .avm_local_bb4_ld__u26_inst0_readdata(avm_local_bb4_ld__u26_inst0_readdata),
      .avm_local_bb4_ld__u26_inst0_readdatavalid(avm_local_bb4_ld__u26_inst0_readdatavalid),
      .avm_local_bb4_ld__u26_inst0_writeack(avm_local_bb4_ld__u26_inst0_writeack),
      // AVM avm_local_bb4_ld__u27_inst0
      .avm_local_bb4_ld__u27_inst0_enable(avm_local_bb4_ld__u27_inst0_enable),
      .avm_local_bb4_ld__u27_inst0_read(avm_local_bb4_ld__u27_inst0_read),
      .avm_local_bb4_ld__u27_inst0_write(avm_local_bb4_ld__u27_inst0_write),
      .avm_local_bb4_ld__u27_inst0_burstcount(avm_local_bb4_ld__u27_inst0_burstcount),
      .avm_local_bb4_ld__u27_inst0_address(avm_local_bb4_ld__u27_inst0_address),
      .avm_local_bb4_ld__u27_inst0_writedata(avm_local_bb4_ld__u27_inst0_writedata),
      .avm_local_bb4_ld__u27_inst0_byteenable(avm_local_bb4_ld__u27_inst0_byteenable),
      .avm_local_bb4_ld__u27_inst0_waitrequest(avm_local_bb4_ld__u27_inst0_waitrequest),
      .avm_local_bb4_ld__u27_inst0_readdata(avm_local_bb4_ld__u27_inst0_readdata),
      .avm_local_bb4_ld__u27_inst0_readdatavalid(avm_local_bb4_ld__u27_inst0_readdatavalid),
      .avm_local_bb4_ld__u27_inst0_writeack(avm_local_bb4_ld__u27_inst0_writeack),
      // AVM avm_local_bb4_ld__u28_inst0
      .avm_local_bb4_ld__u28_inst0_enable(avm_local_bb4_ld__u28_inst0_enable),
      .avm_local_bb4_ld__u28_inst0_read(avm_local_bb4_ld__u28_inst0_read),
      .avm_local_bb4_ld__u28_inst0_write(avm_local_bb4_ld__u28_inst0_write),
      .avm_local_bb4_ld__u28_inst0_burstcount(avm_local_bb4_ld__u28_inst0_burstcount),
      .avm_local_bb4_ld__u28_inst0_address(avm_local_bb4_ld__u28_inst0_address),
      .avm_local_bb4_ld__u28_inst0_writedata(avm_local_bb4_ld__u28_inst0_writedata),
      .avm_local_bb4_ld__u28_inst0_byteenable(avm_local_bb4_ld__u28_inst0_byteenable),
      .avm_local_bb4_ld__u28_inst0_waitrequest(avm_local_bb4_ld__u28_inst0_waitrequest),
      .avm_local_bb4_ld__u28_inst0_readdata(avm_local_bb4_ld__u28_inst0_readdata),
      .avm_local_bb4_ld__u28_inst0_readdatavalid(avm_local_bb4_ld__u28_inst0_readdatavalid),
      .avm_local_bb4_ld__u28_inst0_writeack(avm_local_bb4_ld__u28_inst0_writeack),
      // AVM avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_enable(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_enable),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_read(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_read),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_write(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_write),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_burstcount(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_burstcount),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_address(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_address),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writedata(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writedata),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_byteenable(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_byteenable),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_waitrequest(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_waitrequest),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdata(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdata),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_readdatavalid),
      .avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writeack(avm_local_bb4_ld_memcoalesce_KerConst_load_0_inst0_writeack),
      // AVM avm_local_bb2_st__inst0
      .avm_local_bb2_st__inst0_enable(local_avm_aspace17_enable[0][0]),
      .avm_local_bb2_st__inst0_read(local_avm_aspace17_read[0][0]),
      .avm_local_bb2_st__inst0_write(local_avm_aspace17_write[0][0]),
      .avm_local_bb2_st__inst0_burstcount(local_avm_aspace17_burstcount[0][0]),
      .avm_local_bb2_st__inst0_address(local_avm_aspace17_address[0][0]),
      .avm_local_bb2_st__inst0_writedata(local_avm_aspace17_writedata[0][0]),
      .avm_local_bb2_st__inst0_byteenable(local_avm_aspace17_byteenable[0][0]),
      .avm_local_bb2_st__inst0_waitrequest(local_avm_aspace17_waitrequest[0][0]),
      .avm_local_bb2_st__inst0_readdata(local_avm_aspace17_readdata[0][0]),
      .avm_local_bb2_st__inst0_readdatavalid(local_avm_aspace17_readdatavalid[0][0]),
      .avm_local_bb2_st__inst0_writeack(local_avm_aspace17_writeack[0][0]),
      // AVM avm_local_bb4_ld__inst0
      .avm_local_bb4_ld__inst0_enable(local_avm_aspace17_enable[0][1]),
      .avm_local_bb4_ld__inst0_read(local_avm_aspace17_read[0][1]),
      .avm_local_bb4_ld__inst0_write(local_avm_aspace17_write[0][1]),
      .avm_local_bb4_ld__inst0_burstcount(local_avm_aspace17_burstcount[0][1]),
      .avm_local_bb4_ld__inst0_address(local_avm_aspace17_address[0][1]),
      .avm_local_bb4_ld__inst0_writedata(local_avm_aspace17_writedata[0][1]),
      .avm_local_bb4_ld__inst0_byteenable(local_avm_aspace17_byteenable[0][1]),
      .avm_local_bb4_ld__inst0_waitrequest(local_avm_aspace17_waitrequest[0][1]),
      .avm_local_bb4_ld__inst0_readdata(local_avm_aspace17_readdata[0][1]),
      .avm_local_bb4_ld__inst0_readdatavalid(local_avm_aspace17_readdatavalid[0][1]),
      .avm_local_bb4_ld__inst0_writeack(local_avm_aspace17_writeack[0][1]),
      // AVM avm_local_bb4_ld__u14_inst0
      .avm_local_bb4_ld__u14_inst0_enable(local_avm_aspace17_enable[0][2]),
      .avm_local_bb4_ld__u14_inst0_read(local_avm_aspace17_read[0][2]),
      .avm_local_bb4_ld__u14_inst0_write(local_avm_aspace17_write[0][2]),
      .avm_local_bb4_ld__u14_inst0_burstcount(local_avm_aspace17_burstcount[0][2]),
      .avm_local_bb4_ld__u14_inst0_address(local_avm_aspace17_address[0][2]),
      .avm_local_bb4_ld__u14_inst0_writedata(local_avm_aspace17_writedata[0][2]),
      .avm_local_bb4_ld__u14_inst0_byteenable(local_avm_aspace17_byteenable[0][2]),
      .avm_local_bb4_ld__u14_inst0_waitrequest(local_avm_aspace17_waitrequest[0][2]),
      .avm_local_bb4_ld__u14_inst0_readdata(local_avm_aspace17_readdata[0][2]),
      .avm_local_bb4_ld__u14_inst0_readdatavalid(local_avm_aspace17_readdatavalid[0][2]),
      .avm_local_bb4_ld__u14_inst0_writeack(local_avm_aspace17_writeack[0][2]),
      // AVM avm_local_bb2_st__u3_inst0
      .avm_local_bb2_st__u3_inst0_enable(local_avm_aspace19_enable[0][0]),
      .avm_local_bb2_st__u3_inst0_read(local_avm_aspace19_read[0][0]),
      .avm_local_bb2_st__u3_inst0_write(local_avm_aspace19_write[0][0]),
      .avm_local_bb2_st__u3_inst0_burstcount(local_avm_aspace19_burstcount[0][0]),
      .avm_local_bb2_st__u3_inst0_address(local_avm_aspace19_address[0][0]),
      .avm_local_bb2_st__u3_inst0_writedata(local_avm_aspace19_writedata[0][0]),
      .avm_local_bb2_st__u3_inst0_byteenable(local_avm_aspace19_byteenable[0][0]),
      .avm_local_bb2_st__u3_inst0_waitrequest(local_avm_aspace19_waitrequest[0][0]),
      .avm_local_bb2_st__u3_inst0_readdata(local_avm_aspace19_readdata[0][0]),
      .avm_local_bb2_st__u3_inst0_readdatavalid(local_avm_aspace19_readdatavalid[0][0]),
      .avm_local_bb2_st__u3_inst0_writeack(local_avm_aspace19_writeack[0][0]),
      // AVM avm_local_bb4_ld__u13_inst0
      .avm_local_bb4_ld__u13_inst0_enable(local_avm_aspace19_enable[0][1]),
      .avm_local_bb4_ld__u13_inst0_read(local_avm_aspace19_read[0][1]),
      .avm_local_bb4_ld__u13_inst0_write(local_avm_aspace19_write[0][1]),
      .avm_local_bb4_ld__u13_inst0_burstcount(local_avm_aspace19_burstcount[0][1]),
      .avm_local_bb4_ld__u13_inst0_address(local_avm_aspace19_address[0][1]),
      .avm_local_bb4_ld__u13_inst0_writedata(local_avm_aspace19_writedata[0][1]),
      .avm_local_bb4_ld__u13_inst0_byteenable(local_avm_aspace19_byteenable[0][1]),
      .avm_local_bb4_ld__u13_inst0_waitrequest(local_avm_aspace19_waitrequest[0][1]),
      .avm_local_bb4_ld__u13_inst0_readdata(local_avm_aspace19_readdata[0][1]),
      .avm_local_bb4_ld__u13_inst0_readdatavalid(local_avm_aspace19_readdatavalid[0][1]),
      .avm_local_bb4_ld__u13_inst0_writeack(local_avm_aspace19_writeack[0][1]),
      // AVM avm_local_bb4_ld__u15_inst0
      .avm_local_bb4_ld__u15_inst0_enable(local_avm_aspace19_enable[0][2]),
      .avm_local_bb4_ld__u15_inst0_read(local_avm_aspace19_read[0][2]),
      .avm_local_bb4_ld__u15_inst0_write(local_avm_aspace19_write[0][2]),
      .avm_local_bb4_ld__u15_inst0_burstcount(local_avm_aspace19_burstcount[0][2]),
      .avm_local_bb4_ld__u15_inst0_address(local_avm_aspace19_address[0][2]),
      .avm_local_bb4_ld__u15_inst0_writedata(local_avm_aspace19_writedata[0][2]),
      .avm_local_bb4_ld__u15_inst0_byteenable(local_avm_aspace19_byteenable[0][2]),
      .avm_local_bb4_ld__u15_inst0_waitrequest(local_avm_aspace19_waitrequest[0][2]),
      .avm_local_bb4_ld__u15_inst0_readdata(local_avm_aspace19_readdata[0][2]),
      .avm_local_bb4_ld__u15_inst0_readdatavalid(local_avm_aspace19_readdatavalid[0][2]),
      .avm_local_bb4_ld__u15_inst0_writeack(local_avm_aspace19_writeack[0][2]),
      // AVM avm_local_bb2_st__u5_inst0
      .avm_local_bb2_st__u5_inst0_enable(local_avm_aspace21_enable[0][0]),
      .avm_local_bb2_st__u5_inst0_read(local_avm_aspace21_read[0][0]),
      .avm_local_bb2_st__u5_inst0_write(local_avm_aspace21_write[0][0]),
      .avm_local_bb2_st__u5_inst0_burstcount(local_avm_aspace21_burstcount[0][0]),
      .avm_local_bb2_st__u5_inst0_address(local_avm_aspace21_address[0][0]),
      .avm_local_bb2_st__u5_inst0_writedata(local_avm_aspace21_writedata[0][0]),
      .avm_local_bb2_st__u5_inst0_byteenable(local_avm_aspace21_byteenable[0][0]),
      .avm_local_bb2_st__u5_inst0_waitrequest(local_avm_aspace21_waitrequest[0][0]),
      .avm_local_bb2_st__u5_inst0_readdata(local_avm_aspace21_readdata[0][0]),
      .avm_local_bb2_st__u5_inst0_readdatavalid(local_avm_aspace21_readdatavalid[0][0]),
      .avm_local_bb2_st__u5_inst0_writeack(local_avm_aspace21_writeack[0][0]),
      // AVM avm_local_bb4_ld__u16_inst0
      .avm_local_bb4_ld__u16_inst0_enable(local_avm_aspace21_enable[0][1]),
      .avm_local_bb4_ld__u16_inst0_read(local_avm_aspace21_read[0][1]),
      .avm_local_bb4_ld__u16_inst0_write(local_avm_aspace21_write[0][1]),
      .avm_local_bb4_ld__u16_inst0_burstcount(local_avm_aspace21_burstcount[0][1]),
      .avm_local_bb4_ld__u16_inst0_address(local_avm_aspace21_address[0][1]),
      .avm_local_bb4_ld__u16_inst0_writedata(local_avm_aspace21_writedata[0][1]),
      .avm_local_bb4_ld__u16_inst0_byteenable(local_avm_aspace21_byteenable[0][1]),
      .avm_local_bb4_ld__u16_inst0_waitrequest(local_avm_aspace21_waitrequest[0][1]),
      .avm_local_bb4_ld__u16_inst0_readdata(local_avm_aspace21_readdata[0][1]),
      .avm_local_bb4_ld__u16_inst0_readdatavalid(local_avm_aspace21_readdatavalid[0][1]),
      .avm_local_bb4_ld__u16_inst0_writeack(local_avm_aspace21_writeack[0][1]),
      // AVM avm_local_bb4_ld__u17_inst0
      .avm_local_bb4_ld__u17_inst0_enable(local_avm_aspace21_enable[0][2]),
      .avm_local_bb4_ld__u17_inst0_read(local_avm_aspace21_read[0][2]),
      .avm_local_bb4_ld__u17_inst0_write(local_avm_aspace21_write[0][2]),
      .avm_local_bb4_ld__u17_inst0_burstcount(local_avm_aspace21_burstcount[0][2]),
      .avm_local_bb4_ld__u17_inst0_address(local_avm_aspace21_address[0][2]),
      .avm_local_bb4_ld__u17_inst0_writedata(local_avm_aspace21_writedata[0][2]),
      .avm_local_bb4_ld__u17_inst0_byteenable(local_avm_aspace21_byteenable[0][2]),
      .avm_local_bb4_ld__u17_inst0_waitrequest(local_avm_aspace21_waitrequest[0][2]),
      .avm_local_bb4_ld__u17_inst0_readdata(local_avm_aspace21_readdata[0][2]),
      .avm_local_bb4_ld__u17_inst0_readdatavalid(local_avm_aspace21_readdatavalid[0][2]),
      .avm_local_bb4_ld__u17_inst0_writeack(local_avm_aspace21_writeack[0][2]),
      // AVST avst_local_bb2__chan_Conf2Intrae_active_inst0
      .avst_local_bb2__chan_Conf2Intrae_active_inst0_valid(avm_channel_id_chan_Conf2Intrae_active_read_valid),
      .avst_local_bb2__chan_Conf2Intrae_active_inst0_ready(avm_channel_id_chan_Conf2Intrae_active_read_ready),
      .avst_local_bb2__chan_Conf2Intrae_active_inst0_data(avm_channel_id_chan_Conf2Intrae_active_read_data),
      // AVST avst_local_bb2__chan_Conf2Intrae_cnt_inst0
      .avst_local_bb2__chan_Conf2Intrae_cnt_inst0_valid(avm_channel_id_chan_Conf2Intrae_cnt_read_valid),
      .avst_local_bb2__chan_Conf2Intrae_cnt_inst0_ready(avm_channel_id_chan_Conf2Intrae_cnt_read_ready),
      .avst_local_bb2__chan_Conf2Intrae_cnt_inst0_data(avm_channel_id_chan_Conf2Intrae_cnt_read_data),
      // AVST avst_local_bb2__chan_Conf2Intrae_mode_inst0
      .avst_local_bb2__chan_Conf2Intrae_mode_inst0_valid(avm_channel_id_chan_Conf2Intrae_mode_read_valid),
      .avst_local_bb2__chan_Conf2Intrae_mode_inst0_ready(avm_channel_id_chan_Conf2Intrae_mode_read_ready),
      .avst_local_bb2__chan_Conf2Intrae_mode_inst0_data(avm_channel_id_chan_Conf2Intrae_mode_read_data),
      // AVST avst_local_bb2__chan_Conf2Intrae_x_inst0
      .avst_local_bb2__chan_Conf2Intrae_x_inst0_valid(avm_channel_id_chan_Conf2Intrae_x_read_valid),
      .avst_local_bb2__chan_Conf2Intrae_x_inst0_ready(avm_channel_id_chan_Conf2Intrae_x_read_ready),
      .avst_local_bb2__chan_Conf2Intrae_x_inst0_data(avm_channel_id_chan_Conf2Intrae_x_read_data),
      // AVST avst_local_bb2__chan_Conf2Intrae_y_inst0
      .avst_local_bb2__chan_Conf2Intrae_y_inst0_valid(avm_channel_id_chan_Conf2Intrae_y_read_valid),
      .avst_local_bb2__chan_Conf2Intrae_y_inst0_ready(avm_channel_id_chan_Conf2Intrae_y_read_ready),
      .avst_local_bb2__chan_Conf2Intrae_y_inst0_data(avm_channel_id_chan_Conf2Intrae_y_read_data),
      // AVST avst_local_bb2__chan_Conf2Intrae_z_inst0
      .avst_local_bb2__chan_Conf2Intrae_z_inst0_valid(avm_channel_id_chan_Conf2Intrae_z_read_valid),
      .avst_local_bb2__chan_Conf2Intrae_z_inst0_ready(avm_channel_id_chan_Conf2Intrae_z_read_ready),
      .avst_local_bb2__chan_Conf2Intrae_z_inst0_data(avm_channel_id_chan_Conf2Intrae_z_read_data),
      // AVST avst_local_bb5__chan_Intrae2Store_active_inst0
      .avst_local_bb5__chan_Intrae2Store_active_inst0_valid(avm_channel_id_chan_Intrae2Store_active_write_valid),
      .avst_local_bb5__chan_Intrae2Store_active_inst0_ready(avm_channel_id_chan_Intrae2Store_active_write_ready),
      .avst_local_bb5__chan_Intrae2Store_active_inst0_data(avm_channel_id_chan_Intrae2Store_active_write_data),
      .avst_local_bb5__chan_Intrae2Store_active_inst0_almostfull(avm_channel_id_chan_Intrae2Store_active_write_almostfull),
      // AVST avst_local_bb5__chan_Intrae2Store_cnt_inst0
      .avst_local_bb5__chan_Intrae2Store_cnt_inst0_valid(avm_channel_id_chan_Intrae2Store_cnt_write_valid),
      .avst_local_bb5__chan_Intrae2Store_cnt_inst0_ready(avm_channel_id_chan_Intrae2Store_cnt_write_ready),
      .avst_local_bb5__chan_Intrae2Store_cnt_inst0_data(avm_channel_id_chan_Intrae2Store_cnt_write_data),
      .avst_local_bb5__chan_Intrae2Store_cnt_inst0_almostfull(avm_channel_id_chan_Intrae2Store_cnt_write_almostfull),
      // AVST avst_local_bb5__chan_Intrae2Store_intrae_inst0
      .avst_local_bb5__chan_Intrae2Store_intrae_inst0_valid(avm_channel_id_chan_Intrae2Store_intrae_write_valid),
      .avst_local_bb5__chan_Intrae2Store_intrae_inst0_ready(avm_channel_id_chan_Intrae2Store_intrae_write_ready),
      .avst_local_bb5__chan_Intrae2Store_intrae_inst0_data(avm_channel_id_chan_Intrae2Store_intrae_write_data),
      .avst_local_bb5__chan_Intrae2Store_intrae_inst0_almostfull(avm_channel_id_chan_Intrae2Store_intrae_write_almostfull),
      // AVST avst_local_bb5__chan_Intrae2Store_mode_inst0
      .avst_local_bb5__chan_Intrae2Store_mode_inst0_valid(avm_channel_id_chan_Intrae2Store_mode_write_valid),
      .avst_local_bb5__chan_Intrae2Store_mode_inst0_ready(avm_channel_id_chan_Intrae2Store_mode_write_ready),
      .avst_local_bb5__chan_Intrae2Store_mode_inst0_data(avm_channel_id_chan_Intrae2Store_mode_write_data),
      .avst_local_bb5__chan_Intrae2Store_mode_inst0_almostfull(avm_channel_id_chan_Intrae2Store_mode_write_almostfull),
      // AVM p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb3_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid)
   );

   assign lmem_invalid_single_bit = |lmem_invalid_aspaces;
   generate
   begin:local_mem_system_aspace17
      logic local_icm_arb_request [1][3];
      logic local_icm_arb_enable [1][3];
      logic local_icm_arb_read [1][3];
      logic local_icm_arb_write [1][3];
      logic local_icm_arb_burstcount [1][3];
      logic [6:0] local_icm_arb_address [1][3];
      logic [31:0] local_icm_arb_writedata [1][3];
      logic [3:0] local_icm_arb_byteenable [1][3];
      logic local_icm_arb_stall [1][3];
      logic local_icm_wrp_ack [1][3];
      logic local_icm_rrp_datavalid [1][3];
      logic [31:0] local_icm_rrp_data [1][3];

      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:local_mem_group
         for( __j = 0; __j < 3; __j = __j + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace17_enable[__i][__j]),
               .avm_read(local_avm_aspace17_read[__i][__j]),
               .avm_write(local_avm_aspace17_write[__i][__j]),
               .avm_burstcount(local_avm_aspace17_burstcount[__i][__j]),
               .avm_address(local_avm_aspace17_address[__i][__j]),
               .avm_writedata(local_avm_aspace17_writedata[__i][__j]),
               .avm_byteenable(local_avm_aspace17_byteenable[__i][__j]),
               .avm_waitrequest(local_avm_aspace17_waitrequest[__i][__j]),
               .avm_readdata(local_avm_aspace17_readdata[__i][__j]),
               .avm_readdatavalid(local_avm_aspace17_readdatavalid[__i][__j]),
               .avm_writeack(local_avm_aspace17_writeack[__i][__j]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__i][__j]),
               .ic_arb_enable(local_icm_arb_enable[__i][__j]),
               .ic_arb_read(local_icm_arb_read[__i][__j]),
               .ic_arb_write(local_icm_arb_write[__i][__j]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__i][__j]),
               .ic_arb_address(local_icm_arb_address[__i][__j]),
               .ic_arb_writedata(local_icm_arb_writedata[__i][__j]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__i][__j]),
               .ic_arb_stall(local_icm_arb_stall[__i][__j]),
               .ic_wrp_ack(local_icm_wrp_ack[__i][__j]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__i][__j]),
               .ic_rrp_data(local_icm_rrp_data[__i][__j])
            );

         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:bank
            logic port_enable [1:4];
            logic port_read [1:4];
            logic port_write [1:4];
            logic [6:0] port_address [1:4];
            logic [31:0] port_writedata [1:4];
            logic [3:0] port_byteenable [1:4];
            logic port_waitrequest [1:4];
            logic [31:0] port_readdata [1:4];
            logic port_readdatavalid [1:4];

            // INST mem0 of acl_mem2x
            acl_mem2x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("BIDIR_DUAL_PORT"),
               .PREFERRED_WIDTH(160),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .clk2x(clock2x),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2]),
               // AVS avs_port3
               .avs_port3_enable(port_enable[3]),
               .avs_port3_read(port_read[3]),
               .avs_port3_write(port_write[3]),
               .avs_port3_address(port_address[3]),
               .avs_port3_writedata(port_writedata[3]),
               .avs_port3_byteenable(port_byteenable[3]),
               .avs_port3_waitrequest(port_waitrequest[3]),
               .avs_port3_readdata(port_readdata[3]),
               .avs_port3_readdatavalid(port_readdatavalid[3]),
               // AVS avs_port4
               .avs_port4_enable(port_enable[4]),
               .avs_port4_read(port_read[4]),
               .avs_port4_write(port_write[4]),
               .avs_port4_address(port_address[4]),
               .avs_port4_writedata(port_writedata[4]),
               .avs_port4_byteenable(port_byteenable[4]),
               .avs_port4_waitrequest(port_waitrequest[4]),
               .avs_port4_readdata(port_readdata[4]),
               .avs_port4_readdatavalid(port_readdatavalid[4])
            );

         end

         for( __j = 0; __j < 3; __j = __j + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__i][__j]),
               .m_arb_enable(local_icm_arb_enable[__i][__j]),
               .m_arb_read(local_icm_arb_read[__i][__j]),
               .m_arb_write(local_icm_arb_write[__i][__j]),
               .m_arb_burstcount(local_icm_arb_burstcount[__i][__j]),
               .m_arb_address(local_icm_arb_address[__i][__j]),
               .m_arb_writedata(local_icm_arb_writedata[__i][__j]),
               .m_arb_byteenable(local_icm_arb_byteenable[__i][__j]),
               .m_arb_stall(local_icm_arb_stall[__i][__j]),
               .m_wrp_ack(local_icm_wrp_ack[__i][__j]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__i][__j]),
               .m_rrp_data(local_icm_rrp_data[__i][__j]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_3
            Krnl_GA_system_interconnect_3 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port3bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[2].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[2].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[2].b_arb_read[0];
            assign icm_in_arb_write[0] = router[2].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[2].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[2].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[2].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[2].b_arb_byteenable[0];
            assign router[2].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[2].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[2].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[2].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[3] = icm_out_arb_enable;
            assign bank[0].port_read[3] = icm_out_arb_read;
            assign bank[0].port_write[3] = icm_out_arb_write;
            assign bank[0].port_address[3] = icm_out_arb_address;
            assign bank[0].port_writedata[3] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[3];
            assign icm_out_rrp_data = bank[0].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __j = 0; __j < 1; __j = __j + 1 )
         begin:port4bank0
            assign bank[0].port_enable[4] = '0;
            assign bank[0].port_read[4] = '0;
            assign bank[0].port_write[4] = '0;
            assign bank[0].port_address[4] = '0;
            assign bank[0].port_writedata[4] = '0;
            assign bank[0].port_byteenable[4] = '0;
         end

      end

   end
   endgenerate

   generate
   begin:local_mem_system_aspace19
      logic local_icm_arb_request [1][3];
      logic local_icm_arb_enable [1][3];
      logic local_icm_arb_read [1][3];
      logic local_icm_arb_write [1][3];
      logic local_icm_arb_burstcount [1][3];
      logic [6:0] local_icm_arb_address [1][3];
      logic [31:0] local_icm_arb_writedata [1][3];
      logic [3:0] local_icm_arb_byteenable [1][3];
      logic local_icm_arb_stall [1][3];
      logic local_icm_wrp_ack [1][3];
      logic local_icm_rrp_datavalid [1][3];
      logic [31:0] local_icm_rrp_data [1][3];

      for( __j = 0; __j < 1; __j = __j + 1 )
      begin:local_mem_group
         for( __k = 0; __k < 3; __k = __k + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace19_enable[__j][__k]),
               .avm_read(local_avm_aspace19_read[__j][__k]),
               .avm_write(local_avm_aspace19_write[__j][__k]),
               .avm_burstcount(local_avm_aspace19_burstcount[__j][__k]),
               .avm_address(local_avm_aspace19_address[__j][__k]),
               .avm_writedata(local_avm_aspace19_writedata[__j][__k]),
               .avm_byteenable(local_avm_aspace19_byteenable[__j][__k]),
               .avm_waitrequest(local_avm_aspace19_waitrequest[__j][__k]),
               .avm_readdata(local_avm_aspace19_readdata[__j][__k]),
               .avm_readdatavalid(local_avm_aspace19_readdatavalid[__j][__k]),
               .avm_writeack(local_avm_aspace19_writeack[__j][__k]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__j][__k]),
               .ic_arb_enable(local_icm_arb_enable[__j][__k]),
               .ic_arb_read(local_icm_arb_read[__j][__k]),
               .ic_arb_write(local_icm_arb_write[__j][__k]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__j][__k]),
               .ic_arb_address(local_icm_arb_address[__j][__k]),
               .ic_arb_writedata(local_icm_arb_writedata[__j][__k]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__j][__k]),
               .ic_arb_stall(local_icm_arb_stall[__j][__k]),
               .ic_wrp_ack(local_icm_wrp_ack[__j][__k]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__j][__k]),
               .ic_rrp_data(local_icm_rrp_data[__j][__k])
            );

         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:bank
            logic port_enable [1:4];
            logic port_read [1:4];
            logic port_write [1:4];
            logic [6:0] port_address [1:4];
            logic [31:0] port_writedata [1:4];
            logic [3:0] port_byteenable [1:4];
            logic port_waitrequest [1:4];
            logic [31:0] port_readdata [1:4];
            logic port_readdatavalid [1:4];

            // INST mem0 of acl_mem2x
            acl_mem2x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("BIDIR_DUAL_PORT"),
               .PREFERRED_WIDTH(160),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .clk2x(clock2x),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2]),
               // AVS avs_port3
               .avs_port3_enable(port_enable[3]),
               .avs_port3_read(port_read[3]),
               .avs_port3_write(port_write[3]),
               .avs_port3_address(port_address[3]),
               .avs_port3_writedata(port_writedata[3]),
               .avs_port3_byteenable(port_byteenable[3]),
               .avs_port3_waitrequest(port_waitrequest[3]),
               .avs_port3_readdata(port_readdata[3]),
               .avs_port3_readdatavalid(port_readdatavalid[3]),
               // AVS avs_port4
               .avs_port4_enable(port_enable[4]),
               .avs_port4_read(port_read[4]),
               .avs_port4_write(port_write[4]),
               .avs_port4_address(port_address[4]),
               .avs_port4_writedata(port_writedata[4]),
               .avs_port4_byteenable(port_byteenable[4]),
               .avs_port4_waitrequest(port_waitrequest[4]),
               .avs_port4_readdata(port_readdata[4]),
               .avs_port4_readdatavalid(port_readdatavalid[4])
            );

         end

         for( __k = 0; __k < 3; __k = __k + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__j][__k]),
               .m_arb_enable(local_icm_arb_enable[__j][__k]),
               .m_arb_read(local_icm_arb_read[__j][__k]),
               .m_arb_write(local_icm_arb_write[__j][__k]),
               .m_arb_burstcount(local_icm_arb_burstcount[__j][__k]),
               .m_arb_address(local_icm_arb_address[__j][__k]),
               .m_arb_writedata(local_icm_arb_writedata[__j][__k]),
               .m_arb_byteenable(local_icm_arb_byteenable[__j][__k]),
               .m_arb_stall(local_icm_arb_stall[__j][__k]),
               .m_wrp_ack(local_icm_wrp_ack[__j][__k]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__j][__k]),
               .m_rrp_data(local_icm_rrp_data[__j][__k]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_3
            Krnl_GA_system_interconnect_3 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port3bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[2].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[2].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[2].b_arb_read[0];
            assign icm_in_arb_write[0] = router[2].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[2].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[2].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[2].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[2].b_arb_byteenable[0];
            assign router[2].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[2].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[2].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[2].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[3] = icm_out_arb_enable;
            assign bank[0].port_read[3] = icm_out_arb_read;
            assign bank[0].port_write[3] = icm_out_arb_write;
            assign bank[0].port_address[3] = icm_out_arb_address;
            assign bank[0].port_writedata[3] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[3];
            assign icm_out_rrp_data = bank[0].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __k = 0; __k < 1; __k = __k + 1 )
         begin:port4bank0
            assign bank[0].port_enable[4] = '0;
            assign bank[0].port_read[4] = '0;
            assign bank[0].port_write[4] = '0;
            assign bank[0].port_address[4] = '0;
            assign bank[0].port_writedata[4] = '0;
            assign bank[0].port_byteenable[4] = '0;
         end

      end

   end
   endgenerate

   generate
   begin:local_mem_system_aspace21
      logic local_icm_arb_request [1][3];
      logic local_icm_arb_enable [1][3];
      logic local_icm_arb_read [1][3];
      logic local_icm_arb_write [1][3];
      logic local_icm_arb_burstcount [1][3];
      logic [6:0] local_icm_arb_address [1][3];
      logic [31:0] local_icm_arb_writedata [1][3];
      logic [3:0] local_icm_arb_byteenable [1][3];
      logic local_icm_arb_stall [1][3];
      logic local_icm_wrp_ack [1][3];
      logic local_icm_rrp_datavalid [1][3];
      logic [31:0] local_icm_rrp_data [1][3];

      for( __k = 0; __k < 1; __k = __k + 1 )
      begin:local_mem_group
         for( __l = 0; __l < 3; __l = __l + 1 )
         begin:master
            // INST avm_to_ic of acl_avm_to_ic
            acl_avm_to_ic
            #(
               .DATA_W(32),
               .WRITEDATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(32),
               .BYTEENA_W(4)
            )
            avm_to_ic
            (
               // AVM avm
               .avm_enable(local_avm_aspace21_enable[__k][__l]),
               .avm_read(local_avm_aspace21_read[__k][__l]),
               .avm_write(local_avm_aspace21_write[__k][__l]),
               .avm_burstcount(local_avm_aspace21_burstcount[__k][__l]),
               .avm_address(local_avm_aspace21_address[__k][__l]),
               .avm_writedata(local_avm_aspace21_writedata[__k][__l]),
               .avm_byteenable(local_avm_aspace21_byteenable[__k][__l]),
               .avm_waitrequest(local_avm_aspace21_waitrequest[__k][__l]),
               .avm_readdata(local_avm_aspace21_readdata[__k][__l]),
               .avm_readdatavalid(local_avm_aspace21_readdatavalid[__k][__l]),
               .avm_writeack(local_avm_aspace21_writeack[__k][__l]),
               // ICM ic
               .ic_arb_request(local_icm_arb_request[__k][__l]),
               .ic_arb_enable(local_icm_arb_enable[__k][__l]),
               .ic_arb_read(local_icm_arb_read[__k][__l]),
               .ic_arb_write(local_icm_arb_write[__k][__l]),
               .ic_arb_burstcount(local_icm_arb_burstcount[__k][__l]),
               .ic_arb_address(local_icm_arb_address[__k][__l]),
               .ic_arb_writedata(local_icm_arb_writedata[__k][__l]),
               .ic_arb_byteenable(local_icm_arb_byteenable[__k][__l]),
               .ic_arb_stall(local_icm_arb_stall[__k][__l]),
               .ic_wrp_ack(local_icm_wrp_ack[__k][__l]),
               .ic_rrp_datavalid(local_icm_rrp_datavalid[__k][__l]),
               .ic_rrp_data(local_icm_rrp_data[__k][__l])
            );

         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:bank
            logic port_enable [1:4];
            logic port_read [1:4];
            logic port_write [1:4];
            logic [6:0] port_address [1:4];
            logic [31:0] port_writedata [1:4];
            logic [3:0] port_byteenable [1:4];
            logic port_waitrequest [1:4];
            logic [31:0] port_readdata [1:4];
            logic port_readdatavalid [1:4];

            // INST mem0 of acl_mem2x
            acl_mem2x
            #(
               .INTENDED_DEVICE_FAMILY("Arria 10"),
               .DEPTH_WORDS(128),
               .WIDTH(32),
               .ENABLED(0),
               .RDW_MODE("DONT_CARE"),
               .RAM_OPERATION_MODE("BIDIR_DUAL_PORT"),
               .PREFERRED_WIDTH(160),
               .RAM_BLOCK_TYPE("M20K")
            )
            mem0
            (
               .clk(clock),
               .clk2x(clock2x),
               .resetn(resetn),
               // AVS avs_port1
               .avs_port1_enable(port_enable[1]),
               .avs_port1_read(port_read[1]),
               .avs_port1_write(port_write[1]),
               .avs_port1_address(port_address[1]),
               .avs_port1_writedata(port_writedata[1]),
               .avs_port1_byteenable(port_byteenable[1]),
               .avs_port1_waitrequest(port_waitrequest[1]),
               .avs_port1_readdata(port_readdata[1]),
               .avs_port1_readdatavalid(port_readdatavalid[1]),
               // AVS avs_port2
               .avs_port2_enable(port_enable[2]),
               .avs_port2_read(port_read[2]),
               .avs_port2_write(port_write[2]),
               .avs_port2_address(port_address[2]),
               .avs_port2_writedata(port_writedata[2]),
               .avs_port2_byteenable(port_byteenable[2]),
               .avs_port2_waitrequest(port_waitrequest[2]),
               .avs_port2_readdata(port_readdata[2]),
               .avs_port2_readdatavalid(port_readdatavalid[2]),
               // AVS avs_port3
               .avs_port3_enable(port_enable[3]),
               .avs_port3_read(port_read[3]),
               .avs_port3_write(port_write[3]),
               .avs_port3_address(port_address[3]),
               .avs_port3_writedata(port_writedata[3]),
               .avs_port3_byteenable(port_byteenable[3]),
               .avs_port3_waitrequest(port_waitrequest[3]),
               .avs_port3_readdata(port_readdata[3]),
               .avs_port3_readdatavalid(port_readdatavalid[3]),
               // AVS avs_port4
               .avs_port4_enable(port_enable[4]),
               .avs_port4_read(port_read[4]),
               .avs_port4_write(port_write[4]),
               .avs_port4_address(port_address[4]),
               .avs_port4_writedata(port_writedata[4]),
               .avs_port4_byteenable(port_byteenable[4]),
               .avs_port4_waitrequest(port_waitrequest[4]),
               .avs_port4_readdata(port_readdata[4]),
               .avs_port4_readdatavalid(port_readdatavalid[4])
            );

         end

         for( __l = 0; __l < 3; __l = __l + 1 )
         begin:router
            logic b_arb_request [1];
            logic b_arb_enable [1];
            logic b_arb_read [1];
            logic b_arb_write [1];
            logic b_arb_burstcount [1];
            logic [6:0] b_arb_address [1];
            logic [31:0] b_arb_writedata [1];
            logic [3:0] b_arb_byteenable [1];
            logic b_arb_stall [1];
            logic b_wrp_ack [1];
            logic b_rrp_datavalid [1];
            logic [31:0] b_rrp_data [1];
            logic bank_select;

            // INST router of acl_ic_local_mem_router
            acl_ic_local_mem_router
            #(
               .DATA_W(32),
               .BURSTCOUNT_W(1),
               .ADDRESS_W(7),
               .BYTEENA_W(4),
               .NUM_BANKS(1)
            )
            router
            (
               .clock(clock),
               .resetn(resetn),
               .bank_select(bank_select),
               // ICM m
               .m_arb_request(local_icm_arb_request[__k][__l]),
               .m_arb_enable(local_icm_arb_enable[__k][__l]),
               .m_arb_read(local_icm_arb_read[__k][__l]),
               .m_arb_write(local_icm_arb_write[__k][__l]),
               .m_arb_burstcount(local_icm_arb_burstcount[__k][__l]),
               .m_arb_address(local_icm_arb_address[__k][__l]),
               .m_arb_writedata(local_icm_arb_writedata[__k][__l]),
               .m_arb_byteenable(local_icm_arb_byteenable[__k][__l]),
               .m_arb_stall(local_icm_arb_stall[__k][__l]),
               .m_wrp_ack(local_icm_wrp_ack[__k][__l]),
               .m_rrp_datavalid(local_icm_rrp_datavalid[__k][__l]),
               .m_rrp_data(local_icm_rrp_data[__k][__l]),
               // ICM b
               .b_arb_request(b_arb_request),
               .b_arb_enable(b_arb_enable),
               .b_arb_read(b_arb_read),
               .b_arb_write(b_arb_write),
               .b_arb_burstcount(b_arb_burstcount),
               .b_arb_address(b_arb_address),
               .b_arb_writedata(b_arb_writedata),
               .b_arb_byteenable(b_arb_byteenable),
               .b_arb_stall(b_arb_stall),
               .b_wrp_ack(b_wrp_ack),
               .b_rrp_datavalid(b_rrp_datavalid),
               .b_rrp_data(b_rrp_data)
            );

            assign bank_select = 1'b1;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port1bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[0].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[0].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[0].b_arb_read[0];
            assign icm_in_arb_write[0] = router[0].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[0].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[0].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[0].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[0].b_arb_byteenable[0];
            assign router[0].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[0].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[0].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[0].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_3
            Krnl_GA_system_interconnect_3 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[1] = icm_out_arb_enable;
            assign bank[0].port_read[1] = icm_out_arb_read;
            assign bank[0].port_write[1] = icm_out_arb_write;
            assign bank[0].port_address[1] = icm_out_arb_address;
            assign bank[0].port_writedata[1] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[1] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[1];
            assign icm_out_rrp_data = bank[0].port_readdata[1];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[1];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port2bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[1].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[1].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[1].b_arb_read[0];
            assign icm_in_arb_write[0] = router[1].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[1].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[1].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[1].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[1].b_arb_byteenable[0];
            assign router[1].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[1].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[1].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[1].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[2] = icm_out_arb_enable;
            assign bank[0].port_read[2] = icm_out_arb_read;
            assign bank[0].port_write[2] = icm_out_arb_write;
            assign bank[0].port_address[2] = icm_out_arb_address;
            assign bank[0].port_writedata[2] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[2] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[2];
            assign icm_out_rrp_data = bank[0].port_readdata[2];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[2];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port3bank0
            logic icm_in_arb_request [1];
            logic icm_in_arb_enable [1];
            logic icm_in_arb_read [1];
            logic icm_in_arb_write [1];
            logic icm_in_arb_burstcount [1];
            logic [6:0] icm_in_arb_address [1];
            logic [31:0] icm_in_arb_writedata [1];
            logic [3:0] icm_in_arb_byteenable [1];
            logic icm_in_arb_stall [1];
            logic icm_in_wrp_ack [1];
            logic icm_in_rrp_datavalid [1];
            logic [31:0] icm_in_rrp_data [1];
            logic icm_out_arb_request;
            logic icm_out_arb_enable;
            logic icm_out_arb_read;
            logic icm_out_arb_write;
            logic icm_out_arb_burstcount;
            logic [6:0] icm_out_arb_address;
            logic [31:0] icm_out_arb_writedata;
            logic [3:0] icm_out_arb_byteenable;
            logic icm_out_arb_stall;
            logic icm_out_wrp_ack;
            logic icm_out_rrp_datavalid;
            logic [31:0] icm_out_rrp_data;

            assign icm_in_arb_request[0] = router[2].b_arb_request[0];
            assign icm_in_arb_enable[0] = router[2].b_arb_enable[0];
            assign icm_in_arb_read[0] = router[2].b_arb_read[0];
            assign icm_in_arb_write[0] = router[2].b_arb_write[0];
            assign icm_in_arb_burstcount[0] = router[2].b_arb_burstcount[0];
            assign icm_in_arb_address[0] = router[2].b_arb_address[0];
            assign icm_in_arb_writedata[0] = router[2].b_arb_writedata[0];
            assign icm_in_arb_byteenable[0] = router[2].b_arb_byteenable[0];
            assign router[2].b_arb_stall[0] = icm_in_arb_stall[0];
            assign router[2].b_wrp_ack[0] = icm_in_wrp_ack[0];
            assign router[2].b_rrp_datavalid[0] = icm_in_rrp_datavalid[0];
            assign router[2].b_rrp_data[0] = icm_in_rrp_data[0];
            // INST data_ic of Krnl_GA_system_interconnect_2
            Krnl_GA_system_interconnect_2 data_ic
            (
               .clock(clock),
               .resetn(resetn),
               // ICM m
               .m_arb_request(icm_in_arb_request),
               .m_arb_enable(icm_in_arb_enable),
               .m_arb_read(icm_in_arb_read),
               .m_arb_write(icm_in_arb_write),
               .m_arb_burstcount(icm_in_arb_burstcount),
               .m_arb_address(icm_in_arb_address),
               .m_arb_writedata(icm_in_arb_writedata),
               .m_arb_byteenable(icm_in_arb_byteenable),
               .m_arb_stall(icm_in_arb_stall),
               .m_wrp_ack(icm_in_wrp_ack),
               .m_rrp_datavalid(icm_in_rrp_datavalid),
               .m_rrp_data(icm_in_rrp_data),
               // ICM mout
               .mout_arb_request(icm_out_arb_request),
               .mout_arb_enable(icm_out_arb_enable),
               .mout_arb_read(icm_out_arb_read),
               .mout_arb_write(icm_out_arb_write),
               .mout_arb_burstcount(icm_out_arb_burstcount),
               .mout_arb_address(icm_out_arb_address),
               .mout_arb_writedata(icm_out_arb_writedata),
               .mout_arb_byteenable(icm_out_arb_byteenable),
               .mout_arb_id(),
               .mout_arb_stall(icm_out_arb_stall),
               .mout_wrp_ack(icm_out_wrp_ack),
               .mout_rrp_datavalid(icm_out_rrp_datavalid),
               .mout_rrp_data(icm_out_rrp_data)
            );

            assign bank[0].port_enable[3] = icm_out_arb_enable;
            assign bank[0].port_read[3] = icm_out_arb_read;
            assign bank[0].port_write[3] = icm_out_arb_write;
            assign bank[0].port_address[3] = icm_out_arb_address;
            assign bank[0].port_writedata[3] = icm_out_arb_writedata;
            assign bank[0].port_byteenable[3] = icm_out_arb_byteenable;
            assign icm_out_arb_stall = bank[0].port_waitrequest[3];
            assign icm_out_rrp_data = bank[0].port_readdata[3];
            assign icm_out_rrp_datavalid = bank[0].port_readdatavalid[3];
            assign icm_out_wrp_ack = 'b0;
         end

         for( __l = 0; __l < 1; __l = __l + 1 )
         begin:port4bank0
            assign bank[0].port_enable[4] = '0;
            assign bank[0].port_read[4] = '0;
            assign bank[0].port_write[4] = '0;
            assign bank[0].port_address[4] = '0;
            assign bank[0].port_writedata[4] = '0;
            assign bank[0].port_byteenable[4] = '0;
         end

      end

   end
   endgenerate

endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_Store_top_wrapper_0
/////////////////////////////////////////////////////////////////
module Krnl_Store_top_wrapper_0
(
   input logic start,
   input logic [223:0] kernel_arguments,
   input logic [31:0] work_dim,
   input logic [31:0] global_offset [2:0],
   output logic kernel_valid_out,
   output logic has_a_write_pending,
   output logic has_a_lsu_active,
   input logic [31:0] global_id [2:0],
   input logic [31:0] local_id [2:0],
   input logic [31:0] group_id [2:0],
   input logic [31:0] global_size [2:0],
   input logic [31:0] local_size [2:0],
   input logic [31:0] num_groups [2:0],
   input logic [31:0] workgroup_size,
   output logic kernel_stall_out,
   input logic kernel_valid_in,
   input logic clock,
   input logic resetn,
   input logic clock2x,
   // AVM avm_local_bb1_st__inst0
   output logic avm_local_bb1_st__inst0_enable,
   output logic avm_local_bb1_st__inst0_read,
   output logic avm_local_bb1_st__inst0_write,
   output logic [4:0] avm_local_bb1_st__inst0_burstcount,
   output logic [30:0] avm_local_bb1_st__inst0_address,
   output logic [511:0] avm_local_bb1_st__inst0_writedata,
   output logic [63:0] avm_local_bb1_st__inst0_byteenable,
   input logic avm_local_bb1_st__inst0_waitrequest,
   input logic [511:0] avm_local_bb1_st__inst0_readdata,
   input logic avm_local_bb1_st__inst0_readdatavalid,
   input logic avm_local_bb1_st__inst0_writeack,
   // AVM avm_local_bb1_st__u5_inst0
   output logic avm_local_bb1_st__u5_inst0_enable,
   output logic avm_local_bb1_st__u5_inst0_read,
   output logic avm_local_bb1_st__u5_inst0_write,
   output logic [4:0] avm_local_bb1_st__u5_inst0_burstcount,
   output logic [30:0] avm_local_bb1_st__u5_inst0_address,
   output logic [511:0] avm_local_bb1_st__u5_inst0_writedata,
   output logic [63:0] avm_local_bb1_st__u5_inst0_byteenable,
   input logic avm_local_bb1_st__u5_inst0_waitrequest,
   input logic [511:0] avm_local_bb1_st__u5_inst0_readdata,
   input logic avm_local_bb1_st__u5_inst0_readdatavalid,
   input logic avm_local_bb1_st__u5_inst0_writeack,
   // AVM avm_local_bb1_st__u7_inst0
   output logic avm_local_bb1_st__u7_inst0_enable,
   output logic avm_local_bb1_st__u7_inst0_read,
   output logic avm_local_bb1_st__u7_inst0_write,
   output logic [4:0] avm_local_bb1_st__u7_inst0_burstcount,
   output logic [30:0] avm_local_bb1_st__u7_inst0_address,
   output logic [511:0] avm_local_bb1_st__u7_inst0_writedata,
   output logic [63:0] avm_local_bb1_st__u7_inst0_byteenable,
   input logic avm_local_bb1_st__u7_inst0_waitrequest,
   input logic [511:0] avm_local_bb1_st__u7_inst0_readdata,
   input logic avm_local_bb1_st__u7_inst0_readdatavalid,
   input logic avm_local_bb1_st__u7_inst0_writeack,
   // AVM avm_local_bb1_st__u9_inst0
   output logic avm_local_bb1_st__u9_inst0_enable,
   output logic avm_local_bb1_st__u9_inst0_read,
   output logic avm_local_bb1_st__u9_inst0_write,
   output logic [4:0] avm_local_bb1_st__u9_inst0_burstcount,
   output logic [30:0] avm_local_bb1_st__u9_inst0_address,
   output logic [511:0] avm_local_bb1_st__u9_inst0_writedata,
   output logic [63:0] avm_local_bb1_st__u9_inst0_byteenable,
   input logic avm_local_bb1_st__u9_inst0_waitrequest,
   input logic [511:0] avm_local_bb1_st__u9_inst0_readdata,
   input logic avm_local_bb1_st__u9_inst0_readdatavalid,
   input logic avm_local_bb1_st__u9_inst0_writeack,
   // AVST avm_channel_id_chan_Intere2Store_active_read
   input logic avm_channel_id_chan_Intere2Store_active_read_valid,
   output logic avm_channel_id_chan_Intere2Store_active_read_ready,
   input logic [7:0] avm_channel_id_chan_Intere2Store_active_read_data,
   // AVST avm_channel_id_chan_Intere2Store_cnt_read
   input logic avm_channel_id_chan_Intere2Store_cnt_read_valid,
   output logic avm_channel_id_chan_Intere2Store_cnt_read_ready,
   input logic [31:0] avm_channel_id_chan_Intere2Store_cnt_read_data,
   // AVST avm_channel_id_chan_Intere2Store_intere_read
   input logic avm_channel_id_chan_Intere2Store_intere_read_valid,
   output logic avm_channel_id_chan_Intere2Store_intere_read_ready,
   input logic [31:0] avm_channel_id_chan_Intere2Store_intere_read_data,
   // AVST avm_channel_id_chan_Intere2Store_mode_read
   input logic avm_channel_id_chan_Intere2Store_mode_read_valid,
   output logic avm_channel_id_chan_Intere2Store_mode_read_ready,
   input logic [7:0] avm_channel_id_chan_Intere2Store_mode_read_data,
   // AVST avm_channel_id_chan_Intrae2Store_active_read
   input logic avm_channel_id_chan_Intrae2Store_active_read_valid,
   output logic avm_channel_id_chan_Intrae2Store_active_read_ready,
   input logic [7:0] avm_channel_id_chan_Intrae2Store_active_read_data,
   // AVST avm_channel_id_chan_Intrae2Store_cnt_read
   input logic avm_channel_id_chan_Intrae2Store_cnt_read_valid,
   output logic avm_channel_id_chan_Intrae2Store_cnt_read_ready,
   input logic [31:0] avm_channel_id_chan_Intrae2Store_cnt_read_data,
   // AVST avm_channel_id_chan_Intrae2Store_intrae_read
   input logic avm_channel_id_chan_Intrae2Store_intrae_read_valid,
   output logic avm_channel_id_chan_Intrae2Store_intrae_read_ready,
   input logic [31:0] avm_channel_id_chan_Intrae2Store_intrae_read_data,
   // AVST avm_channel_id_chan_Intrae2Store_mode_read
   input logic avm_channel_id_chan_Intrae2Store_mode_read_valid,
   output logic avm_channel_id_chan_Intrae2Store_mode_read_ready,
   input logic [7:0] avm_channel_id_chan_Intrae2Store_mode_read_data,
   // AVM p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid,
   // AVM p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdatavalid,
   // AVM p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdatavalid,
   // AVM p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0
   output logic p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_enable,
   output logic p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_read,
   output logic p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_write,
   output logic [5:0] p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_burstcount,
   output logic [31:0] p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_address,
   output logic [255:0] p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_writedata,
   output logic [31:0] p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_byteenable,
   input logic p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_waitrequest,
   input logic [255:0] p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdata,
   input logic p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdatavalid
);
   logic lmem_invalid_single_bit;

   // INST kernel of Krnl_Store_function_wrapper
   Krnl_Store_function_wrapper kernel
   (
      .local_router_hang(lmem_invalid_single_bit),
      .start(start),
      .kernel_arguments(kernel_arguments),
      .work_dim(work_dim),
      .global_offset_0(global_offset[0]),
      .global_offset_1(global_offset[1]),
      .global_offset_2(global_offset[2]),
      .kernel_valid_out(kernel_valid_out),
      .has_a_write_pending(has_a_write_pending),
      .has_a_lsu_active(has_a_lsu_active),
      .global_id_0(global_id[0]),
      .global_id_1(global_id[1]),
      .global_id_2(global_id[2]),
      .local_id_0(local_id[0]),
      .local_id_1(local_id[1]),
      .local_id_2(local_id[2]),
      .group_id_0(group_id[0]),
      .group_id_1(group_id[1]),
      .group_id_2(group_id[2]),
      .global_size_0(global_size[0]),
      .global_size_1(global_size[1]),
      .global_size_2(global_size[2]),
      .local_size_0(local_size[0]),
      .local_size_1(local_size[1]),
      .local_size_2(local_size[2]),
      .num_groups_0(num_groups[0]),
      .num_groups_1(num_groups[1]),
      .num_groups_2(num_groups[2]),
      .workgroup_size(workgroup_size),
      .kernel_stall_out(kernel_stall_out),
      .kernel_valid_in(kernel_valid_in),
      .clock(clock),
      .resetn(resetn),
      .clock2x(clock2x),
      // AVM avm_local_bb1_st__inst0
      .avm_local_bb1_st__inst0_enable(avm_local_bb1_st__inst0_enable),
      .avm_local_bb1_st__inst0_read(avm_local_bb1_st__inst0_read),
      .avm_local_bb1_st__inst0_write(avm_local_bb1_st__inst0_write),
      .avm_local_bb1_st__inst0_burstcount(avm_local_bb1_st__inst0_burstcount),
      .avm_local_bb1_st__inst0_address(avm_local_bb1_st__inst0_address),
      .avm_local_bb1_st__inst0_writedata(avm_local_bb1_st__inst0_writedata),
      .avm_local_bb1_st__inst0_byteenable(avm_local_bb1_st__inst0_byteenable),
      .avm_local_bb1_st__inst0_waitrequest(avm_local_bb1_st__inst0_waitrequest),
      .avm_local_bb1_st__inst0_readdata(avm_local_bb1_st__inst0_readdata),
      .avm_local_bb1_st__inst0_readdatavalid(avm_local_bb1_st__inst0_readdatavalid),
      .avm_local_bb1_st__inst0_writeack(avm_local_bb1_st__inst0_writeack),
      // AVM avm_local_bb1_st__u5_inst0
      .avm_local_bb1_st__u5_inst0_enable(avm_local_bb1_st__u5_inst0_enable),
      .avm_local_bb1_st__u5_inst0_read(avm_local_bb1_st__u5_inst0_read),
      .avm_local_bb1_st__u5_inst0_write(avm_local_bb1_st__u5_inst0_write),
      .avm_local_bb1_st__u5_inst0_burstcount(avm_local_bb1_st__u5_inst0_burstcount),
      .avm_local_bb1_st__u5_inst0_address(avm_local_bb1_st__u5_inst0_address),
      .avm_local_bb1_st__u5_inst0_writedata(avm_local_bb1_st__u5_inst0_writedata),
      .avm_local_bb1_st__u5_inst0_byteenable(avm_local_bb1_st__u5_inst0_byteenable),
      .avm_local_bb1_st__u5_inst0_waitrequest(avm_local_bb1_st__u5_inst0_waitrequest),
      .avm_local_bb1_st__u5_inst0_readdata(avm_local_bb1_st__u5_inst0_readdata),
      .avm_local_bb1_st__u5_inst0_readdatavalid(avm_local_bb1_st__u5_inst0_readdatavalid),
      .avm_local_bb1_st__u5_inst0_writeack(avm_local_bb1_st__u5_inst0_writeack),
      // AVM avm_local_bb1_st__u7_inst0
      .avm_local_bb1_st__u7_inst0_enable(avm_local_bb1_st__u7_inst0_enable),
      .avm_local_bb1_st__u7_inst0_read(avm_local_bb1_st__u7_inst0_read),
      .avm_local_bb1_st__u7_inst0_write(avm_local_bb1_st__u7_inst0_write),
      .avm_local_bb1_st__u7_inst0_burstcount(avm_local_bb1_st__u7_inst0_burstcount),
      .avm_local_bb1_st__u7_inst0_address(avm_local_bb1_st__u7_inst0_address),
      .avm_local_bb1_st__u7_inst0_writedata(avm_local_bb1_st__u7_inst0_writedata),
      .avm_local_bb1_st__u7_inst0_byteenable(avm_local_bb1_st__u7_inst0_byteenable),
      .avm_local_bb1_st__u7_inst0_waitrequest(avm_local_bb1_st__u7_inst0_waitrequest),
      .avm_local_bb1_st__u7_inst0_readdata(avm_local_bb1_st__u7_inst0_readdata),
      .avm_local_bb1_st__u7_inst0_readdatavalid(avm_local_bb1_st__u7_inst0_readdatavalid),
      .avm_local_bb1_st__u7_inst0_writeack(avm_local_bb1_st__u7_inst0_writeack),
      // AVM avm_local_bb1_st__u9_inst0
      .avm_local_bb1_st__u9_inst0_enable(avm_local_bb1_st__u9_inst0_enable),
      .avm_local_bb1_st__u9_inst0_read(avm_local_bb1_st__u9_inst0_read),
      .avm_local_bb1_st__u9_inst0_write(avm_local_bb1_st__u9_inst0_write),
      .avm_local_bb1_st__u9_inst0_burstcount(avm_local_bb1_st__u9_inst0_burstcount),
      .avm_local_bb1_st__u9_inst0_address(avm_local_bb1_st__u9_inst0_address),
      .avm_local_bb1_st__u9_inst0_writedata(avm_local_bb1_st__u9_inst0_writedata),
      .avm_local_bb1_st__u9_inst0_byteenable(avm_local_bb1_st__u9_inst0_byteenable),
      .avm_local_bb1_st__u9_inst0_waitrequest(avm_local_bb1_st__u9_inst0_waitrequest),
      .avm_local_bb1_st__u9_inst0_readdata(avm_local_bb1_st__u9_inst0_readdata),
      .avm_local_bb1_st__u9_inst0_readdatavalid(avm_local_bb1_st__u9_inst0_readdatavalid),
      .avm_local_bb1_st__u9_inst0_writeack(avm_local_bb1_st__u9_inst0_writeack),
      // AVST avst_local_bb1__chan_Intere2Store_active_inst0
      .avst_local_bb1__chan_Intere2Store_active_inst0_valid(avm_channel_id_chan_Intere2Store_active_read_valid),
      .avst_local_bb1__chan_Intere2Store_active_inst0_ready(avm_channel_id_chan_Intere2Store_active_read_ready),
      .avst_local_bb1__chan_Intere2Store_active_inst0_data(avm_channel_id_chan_Intere2Store_active_read_data),
      // AVST avst_local_bb1__chan_Intere2Store_cnt_inst0
      .avst_local_bb1__chan_Intere2Store_cnt_inst0_valid(avm_channel_id_chan_Intere2Store_cnt_read_valid),
      .avst_local_bb1__chan_Intere2Store_cnt_inst0_ready(avm_channel_id_chan_Intere2Store_cnt_read_ready),
      .avst_local_bb1__chan_Intere2Store_cnt_inst0_data(avm_channel_id_chan_Intere2Store_cnt_read_data),
      // AVST avst_local_bb1__chan_Intere2Store_intere_inst0
      .avst_local_bb1__chan_Intere2Store_intere_inst0_valid(avm_channel_id_chan_Intere2Store_intere_read_valid),
      .avst_local_bb1__chan_Intere2Store_intere_inst0_ready(avm_channel_id_chan_Intere2Store_intere_read_ready),
      .avst_local_bb1__chan_Intere2Store_intere_inst0_data(avm_channel_id_chan_Intere2Store_intere_read_data),
      // AVST avst_local_bb1__chan_Intere2Store_mode_inst0
      .avst_local_bb1__chan_Intere2Store_mode_inst0_valid(avm_channel_id_chan_Intere2Store_mode_read_valid),
      .avst_local_bb1__chan_Intere2Store_mode_inst0_ready(avm_channel_id_chan_Intere2Store_mode_read_ready),
      .avst_local_bb1__chan_Intere2Store_mode_inst0_data(avm_channel_id_chan_Intere2Store_mode_read_data),
      // AVST avst_local_bb1__chan_Intrae2Store_active_inst0
      .avst_local_bb1__chan_Intrae2Store_active_inst0_valid(avm_channel_id_chan_Intrae2Store_active_read_valid),
      .avst_local_bb1__chan_Intrae2Store_active_inst0_ready(avm_channel_id_chan_Intrae2Store_active_read_ready),
      .avst_local_bb1__chan_Intrae2Store_active_inst0_data(avm_channel_id_chan_Intrae2Store_active_read_data),
      // AVST avst_local_bb1__chan_Intrae2Store_cnt_inst0
      .avst_local_bb1__chan_Intrae2Store_cnt_inst0_valid(avm_channel_id_chan_Intrae2Store_cnt_read_valid),
      .avst_local_bb1__chan_Intrae2Store_cnt_inst0_ready(avm_channel_id_chan_Intrae2Store_cnt_read_ready),
      .avst_local_bb1__chan_Intrae2Store_cnt_inst0_data(avm_channel_id_chan_Intrae2Store_cnt_read_data),
      // AVST avst_local_bb1__chan_Intrae2Store_intrae_inst0
      .avst_local_bb1__chan_Intrae2Store_intrae_inst0_valid(avm_channel_id_chan_Intrae2Store_intrae_read_valid),
      .avst_local_bb1__chan_Intrae2Store_intrae_inst0_ready(avm_channel_id_chan_Intrae2Store_intrae_read_ready),
      .avst_local_bb1__chan_Intrae2Store_intrae_inst0_data(avm_channel_id_chan_Intrae2Store_intrae_read_data),
      // AVST avst_local_bb1__chan_Intrae2Store_mode_inst0
      .avst_local_bb1__chan_Intrae2Store_mode_inst0_valid(avm_channel_id_chan_Intrae2Store_mode_read_valid),
      .avst_local_bb1__chan_Intrae2Store_mode_inst0_ready(avm_channel_id_chan_Intrae2Store_mode_read_ready),
      .avst_local_bb1__chan_Intrae2Store_mode_inst0_data(avm_channel_id_chan_Intrae2Store_mode_read_data),
      // AVM p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_read(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_write(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_address(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb1_printf_addr_acl_printf_p1i8_32_inst0_readdatavalid),
      // AVM p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_read(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_write(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_address(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb1_printf_addr8_acl_printf_p1i8_32_inst0_readdatavalid),
      // AVM p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_read(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_write(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_address(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb1_printf_addr2_acl_printf_p1i8_32_inst0_readdatavalid),
      // AVM p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_enable(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_enable),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_read(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_read),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_write(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_write),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_burstcount(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_burstcount),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_address(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_address),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_writedata(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_writedata),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_byteenable(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_byteenable),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_waitrequest(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_waitrequest),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdata(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdata),
      .p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdatavalid(p_avm_local_bb1_printf_addr5_acl_printf_p1i8_32_inst0_readdatavalid)
   );

   assign lmem_invalid_single_bit = 'b0;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_0
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_0
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [1],
   input logic m_arb_enable [1],
   input logic m_arb_read [1],
   input logic m_arb_write [1],
   input logic m_arb_burstcount [1],
   input logic [3:0] m_arb_address [1],
   input logic [63:0] m_arb_writedata [1],
   input logic [7:0] m_arb_byteenable [1],
   output logic m_arb_stall [1],
   output logic m_wrp_ack [1],
   output logic m_rrp_datavalid [1],
   output logic [63:0] m_rrp_data [1],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic mout_arb_burstcount,
   output logic [3:0] mout_arb_address,
   output logic [63:0] mout_arb_writedata,
   output logic [7:0] mout_arb_byteenable,
   output logic mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [63:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:m
         logic id;
         acl_ic_master_intf
         #(
            .DATA_W(64),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(4),
            .BYTEENA_W(8),
            .ID_W(1)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(64),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(4),
            .BYTEENA_W(8),
            .ID_W(1)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(1)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(64),
            .ID_W(1)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(64),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(4),
            .BYTEENA_W(8),
            .ID_W(1),
            .TOTAL_NUM_MASTERS(1),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(64),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(4),
         .BYTEENA_W(8),
         .ID_W(1)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(64),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(4),
         .BYTEENA_W(8),
         .ID_W(1)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(1)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(64),
         .ID_W(1)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(64),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(4),
         .BYTEENA_W(8),
         .ID_W(1),
         .NUM_MASTERS(1),
         .PIPELINE_RETURN_PATHS(0),
         .WRP_FIFO_DEPTH(0),
         .RRP_FIFO_DEPTH(0),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(4),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
      assign m[0].wrp_intf.ack = s.wrp_intf.ack;
      assign m[0].wrp_intf.id = s.wrp_intf.id;
   end
   endgenerate

   generate
   begin:rrp
   end
   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = s.in_arb_intf.stall;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_1
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_1
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [1],
   input logic m_arb_enable [1],
   input logic m_arb_read [1],
   input logic m_arb_write [1],
   input logic m_arb_burstcount [1],
   input logic [3:0] m_arb_address [1],
   input logic [63:0] m_arb_writedata [1],
   input logic [7:0] m_arb_byteenable [1],
   output logic m_arb_stall [1],
   output logic m_wrp_ack [1],
   output logic m_rrp_datavalid [1],
   output logic [63:0] m_rrp_data [1],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic mout_arb_burstcount,
   output logic [3:0] mout_arb_address,
   output logic [63:0] mout_arb_writedata,
   output logic [7:0] mout_arb_byteenable,
   output logic mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [63:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:m
         logic id;
         acl_ic_master_intf
         #(
            .DATA_W(64),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(4),
            .BYTEENA_W(8),
            .ID_W(1)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(64),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(4),
            .BYTEENA_W(8),
            .ID_W(1)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(1)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(64),
            .ID_W(1)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(64),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(4),
            .BYTEENA_W(8),
            .ID_W(1),
            .TOTAL_NUM_MASTERS(1),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(64),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(4),
         .BYTEENA_W(8),
         .ID_W(1)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(64),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(4),
         .BYTEENA_W(8),
         .ID_W(1)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(1)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(64),
         .ID_W(1)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(64),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(4),
         .BYTEENA_W(8),
         .ID_W(1),
         .NUM_MASTERS(1),
         .PIPELINE_RETURN_PATHS(0),
         .WRP_FIFO_DEPTH(0),
         .RRP_FIFO_DEPTH(0),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(4),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
   end
   endgenerate

   generate
   begin:rrp
      assign m[0].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[0].rrp_intf.data = s.rrp_intf.data;
      assign m[0].rrp_intf.id = s.rrp_intf.id;
   end
   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = s.in_arb_intf.stall;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_2
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_2
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [1],
   input logic m_arb_enable [1],
   input logic m_arb_read [1],
   input logic m_arb_write [1],
   input logic m_arb_burstcount [1],
   input logic [6:0] m_arb_address [1],
   input logic [31:0] m_arb_writedata [1],
   input logic [3:0] m_arb_byteenable [1],
   output logic m_arb_stall [1],
   output logic m_wrp_ack [1],
   output logic m_rrp_datavalid [1],
   output logic [31:0] m_rrp_data [1],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic mout_arb_burstcount,
   output logic [6:0] mout_arb_address,
   output logic [31:0] mout_arb_writedata,
   output logic [3:0] mout_arb_byteenable,
   output logic mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [31:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:m
         logic id;
         acl_ic_master_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(1)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(32),
            .ID_W(1)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1),
            .TOTAL_NUM_MASTERS(1),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(1)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(32),
         .ID_W(1)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1),
         .NUM_MASTERS(1),
         .PIPELINE_RETURN_PATHS(0),
         .WRP_FIFO_DEPTH(0),
         .RRP_FIFO_DEPTH(0),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(4),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
   end
   endgenerate

   generate
   begin:rrp
      assign m[0].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[0].rrp_intf.data = s.rrp_intf.data;
      assign m[0].rrp_intf.id = s.rrp_intf.id;
   end
   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = s.in_arb_intf.stall;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_3
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_3
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [1],
   input logic m_arb_enable [1],
   input logic m_arb_read [1],
   input logic m_arb_write [1],
   input logic m_arb_burstcount [1],
   input logic [6:0] m_arb_address [1],
   input logic [31:0] m_arb_writedata [1],
   input logic [3:0] m_arb_byteenable [1],
   output logic m_arb_stall [1],
   output logic m_wrp_ack [1],
   output logic m_rrp_datavalid [1],
   output logic [31:0] m_rrp_data [1],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic mout_arb_burstcount,
   output logic [6:0] mout_arb_address,
   output logic [31:0] mout_arb_writedata,
   output logic [3:0] mout_arb_byteenable,
   output logic mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [31:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:m
         logic id;
         acl_ic_master_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(1)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(32),
            .ID_W(1)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1),
            .TOTAL_NUM_MASTERS(1),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(1)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(32),
         .ID_W(1)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1),
         .NUM_MASTERS(1),
         .PIPELINE_RETURN_PATHS(0),
         .WRP_FIFO_DEPTH(0),
         .RRP_FIFO_DEPTH(0),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(4),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
      assign m[0].wrp_intf.ack = s.wrp_intf.ack;
      assign m[0].wrp_intf.id = s.wrp_intf.id;
   end
   endgenerate

   generate
   begin:rrp
   end
   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = s.in_arb_intf.stall;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_4
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_4
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [1],
   input logic m_arb_enable [1],
   input logic m_arb_read [1],
   input logic m_arb_write [1],
   input logic m_arb_burstcount [1],
   input logic [6:0] m_arb_address [1],
   input logic [31:0] m_arb_writedata [1],
   input logic [3:0] m_arb_byteenable [1],
   output logic m_arb_stall [1],
   output logic m_wrp_ack [1],
   output logic m_rrp_datavalid [1],
   output logic [31:0] m_rrp_data [1],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic mout_arb_burstcount,
   output logic [6:0] mout_arb_address,
   output logic [31:0] mout_arb_writedata,
   output logic [3:0] mout_arb_byteenable,
   output logic mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [31:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:m
         logic id;
         acl_ic_master_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(1)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(32),
            .ID_W(1)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1),
            .TOTAL_NUM_MASTERS(1),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(1)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(32),
         .ID_W(1)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1),
         .NUM_MASTERS(1),
         .PIPELINE_RETURN_PATHS(0),
         .WRP_FIFO_DEPTH(0),
         .RRP_FIFO_DEPTH(0),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(3),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
      assign m[0].wrp_intf.ack = s.wrp_intf.ack;
      assign m[0].wrp_intf.id = s.wrp_intf.id;
   end
   endgenerate

   generate
   begin:rrp
   end
   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = s.in_arb_intf.stall;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_5
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_5
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [1],
   input logic m_arb_enable [1],
   input logic m_arb_read [1],
   input logic m_arb_write [1],
   input logic m_arb_burstcount [1],
   input logic [6:0] m_arb_address [1],
   input logic [31:0] m_arb_writedata [1],
   input logic [3:0] m_arb_byteenable [1],
   output logic m_arb_stall [1],
   output logic m_wrp_ack [1],
   output logic m_rrp_datavalid [1],
   output logic [31:0] m_rrp_data [1],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic mout_arb_burstcount,
   output logic [6:0] mout_arb_address,
   output logic [31:0] mout_arb_writedata,
   output logic [3:0] mout_arb_byteenable,
   output logic mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [31:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:m
         logic id;
         acl_ic_master_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(1)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(32),
            .ID_W(1)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(32),
            .BURSTCOUNT_W(1),
            .ADDRESS_W(7),
            .BYTEENA_W(4),
            .ID_W(1),
            .TOTAL_NUM_MASTERS(1),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(1)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(32),
         .ID_W(1)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(32),
         .BURSTCOUNT_W(1),
         .ADDRESS_W(7),
         .BYTEENA_W(4),
         .ID_W(1),
         .NUM_MASTERS(1),
         .PIPELINE_RETURN_PATHS(0),
         .WRP_FIFO_DEPTH(0),
         .RRP_FIFO_DEPTH(0),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(3),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
   end
   endgenerate

   generate
   begin:rrp
      assign m[0].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[0].rrp_intf.data = s.rrp_intf.data;
      assign m[0].rrp_intf.id = s.rrp_intf.id;
   end
   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = s.in_arb_intf.stall;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_6
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_6
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [1],
   input logic m_arb_enable [1],
   input logic m_arb_read [1],
   input logic m_arb_write [1],
   input logic [5:0] m_arb_burstcount [1],
   input logic [26:0] m_arb_address [1],
   input logic [255:0] m_arb_writedata [1],
   input logic [31:0] m_arb_byteenable [1],
   output logic m_arb_stall [1],
   output logic m_wrp_ack [1],
   output logic m_rrp_datavalid [1],
   output logic [255:0] m_rrp_data [1],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic [5:0] mout_arb_burstcount,
   output logic [26:0] mout_arb_address,
   output logic [255:0] mout_arb_writedata,
   output logic [31:0] mout_arb_byteenable,
   output logic mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [255:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:m
         logic id;
         acl_ic_master_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(1)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(1)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(1)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(256),
            .ID_W(1)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(1),
            .TOTAL_NUM_MASTERS(1),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(256),
         .BURSTCOUNT_W(6),
         .ADDRESS_W(27),
         .BYTEENA_W(32),
         .ID_W(1)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(256),
         .BURSTCOUNT_W(6),
         .ADDRESS_W(27),
         .BYTEENA_W(32),
         .ID_W(1)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(1)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(256),
         .ID_W(1)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(256),
         .BURSTCOUNT_W(6),
         .ADDRESS_W(27),
         .BYTEENA_W(32),
         .ID_W(1),
         .NUM_MASTERS(1),
         .PIPELINE_RETURN_PATHS(1),
         .WRP_FIFO_DEPTH(64),
         .RRP_FIFO_DEPTH(64),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(0),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
   end
   endgenerate

   generate
   begin:rrp
      assign m[0].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[0].rrp_intf.data = s.rrp_intf.data;
      assign m[0].rrp_intf.id = s.rrp_intf.id;
   end
   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = s.in_arb_intf.stall;
endmodule

/////////////////////////////////////////////////////////////////
// MODULE Krnl_GA_system_interconnect_7
/////////////////////////////////////////////////////////////////
module Krnl_GA_system_interconnect_7
(
   input logic clock,
   input logic resetn,
   // ICM m
   input logic m_arb_request [4],
   input logic m_arb_enable [4],
   input logic m_arb_read [4],
   input logic m_arb_write [4],
   input logic [5:0] m_arb_burstcount [4],
   input logic [26:0] m_arb_address [4],
   input logic [255:0] m_arb_writedata [4],
   input logic [31:0] m_arb_byteenable [4],
   output logic m_arb_stall [4],
   output logic m_wrp_ack [4],
   output logic m_rrp_datavalid [4],
   output logic [255:0] m_rrp_data [4],
   // ICM mout
   output logic mout_arb_request,
   output logic mout_arb_enable,
   output logic mout_arb_read,
   output logic mout_arb_write,
   output logic [5:0] mout_arb_burstcount,
   output logic [26:0] mout_arb_address,
   output logic [255:0] mout_arb_writedata,
   output logic [31:0] mout_arb_byteenable,
   output logic [1:0] mout_arb_id,
   input logic mout_arb_stall,
   input logic mout_wrp_ack,
   input logic mout_rrp_datavalid,
   input logic [255:0] mout_rrp_data
);
   genvar __i;
   generate
      for( __i = 0; __i < 4; __i = __i + 1 )
      begin:m
         logic [1:0] id;
         acl_ic_master_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         ) m_intf();
         acl_arb_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         ) arb_intf();
         acl_ic_wrp_intf
         #(
            .ID_W(2)
         ) wrp_intf();
         acl_ic_rrp_intf
         #(
            .DATA_W(256),
            .ID_W(2)
         ) rrp_intf();

         assign id = __i;
         // INST m_endp of acl_ic_master_endpoint
         acl_ic_master_endpoint
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2),
            .TOTAL_NUM_MASTERS(4),
            .ID(__i)
         )
         m_endp
         (
            .clock(clock),
            .resetn(resetn),
            .m_intf(m_intf),
            .arb_intf(arb_intf),
            .wrp_intf(wrp_intf),
            .rrp_intf(rrp_intf)
         );

         assign m_intf.arb.req.request = m_arb_request[__i];
         assign m_intf.arb.req.enable = m_arb_enable[__i];
         assign m_intf.arb.req.read = m_arb_read[__i];
         assign m_intf.arb.req.write = m_arb_write[__i];
         assign m_intf.arb.req.burstcount = m_arb_burstcount[__i];
         assign m_intf.arb.req.address = m_arb_address[__i];
         assign m_intf.arb.req.writedata = m_arb_writedata[__i];
         assign m_intf.arb.req.byteenable = m_arb_byteenable[__i];
         assign m_arb_stall[__i] = m_intf.arb.stall;
         assign m_wrp_ack[__i] = m_intf.wrp.ack;
         assign m_rrp_datavalid[__i] = m_intf.rrp.datavalid;
         assign m_rrp_data[__i] = m_intf.rrp.data;
         assign m_intf.arb.req.id = id;
      end

   endgenerate

   generate
   begin:s
      acl_arb_intf
      #(
         .DATA_W(256),
         .BURSTCOUNT_W(6),
         .ADDRESS_W(27),
         .BYTEENA_W(32),
         .ID_W(2)
      ) in_arb_intf();
      acl_arb_intf
      #(
         .DATA_W(256),
         .BURSTCOUNT_W(6),
         .ADDRESS_W(27),
         .BYTEENA_W(32),
         .ID_W(2)
      ) out_arb_intf();
      acl_ic_wrp_intf
      #(
         .ID_W(2)
      ) wrp_intf();
      acl_ic_rrp_intf
      #(
         .DATA_W(256),
         .ID_W(2)
      ) rrp_intf();

      // INST s_endp of acl_ic_slave_endpoint
      acl_ic_slave_endpoint
      #(
         .DATA_W(256),
         .BURSTCOUNT_W(6),
         .ADDRESS_W(27),
         .BYTEENA_W(32),
         .ID_W(2),
         .NUM_MASTERS(4),
         .PIPELINE_RETURN_PATHS(1),
         .WRP_FIFO_DEPTH(64),
         .RRP_FIFO_DEPTH(64),
         .RRP_USE_LL_FIFO(1),
         .SLAVE_FIXED_LATENCY(0),
         .SEPARATE_READ_WRITE_STALLS(0)
      )
      s_endp
      (
         .clock(clock),
         .resetn(resetn),
         .m_intf(in_arb_intf),
         .s_intf(out_arb_intf),
         .s_readdatavalid(mout_rrp_datavalid),
         .s_readdata(mout_rrp_data),
         .s_writeack(mout_wrp_ack),
         .wrp_intf(wrp_intf),
         .rrp_intf(rrp_intf)
      );

   end
   endgenerate

   generate
   begin:wrp
   end
   endgenerate

   generate
   begin:rrp
      assign m[0].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[0].rrp_intf.data = s.rrp_intf.data;
      assign m[0].rrp_intf.id = s.rrp_intf.id;
      assign m[1].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[1].rrp_intf.data = s.rrp_intf.data;
      assign m[1].rrp_intf.id = s.rrp_intf.id;
      assign m[2].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[2].rrp_intf.data = s.rrp_intf.data;
      assign m[2].rrp_intf.id = s.rrp_intf.id;
      assign m[3].rrp_intf.datavalid = s.rrp_intf.datavalid;
      assign m[3].rrp_intf.data = s.rrp_intf.data;
      assign m[3].rrp_intf.id = s.rrp_intf.id;
   end
   endgenerate

   generate
      for( __i = 0; __i < 3; __i = __i + 1 )
      begin:a
         acl_arb_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         ) m0_intf();
         acl_arb_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         ) m1_intf();
         acl_arb_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         ) mout_intf();

         // INST a of acl_arb2
         acl_arb2
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2),
            .PIPELINE("none"),
            .KEEP_LAST_GRANT(0),
            .NO_STALL_NETWORK(0)
         )
         a
         (
            .clock(clock),
            .resetn(resetn),
            .m0_intf(m0_intf),
            .m1_intf(m1_intf),
            .mout_intf(mout_intf)
         );

      end

   endgenerate

   generate
      for( __i = 0; __i < 1; __i = __i + 1 )
      begin:dp
         acl_arb_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         ) in_intf();
         acl_arb_intf
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         ) out_intf();

         // INST dp of acl_arb_pipeline_reg
         acl_arb_pipeline_reg
         #(
            .DATA_W(256),
            .BURSTCOUNT_W(6),
            .ADDRESS_W(27),
            .BYTEENA_W(32),
            .ID_W(2)
         )
         dp
         (
            .clock(clock),
            .resetn(resetn),
            .in_intf(in_intf),
            .out_intf(out_intf)
         );

      end

   endgenerate

   assign mout_arb_request = s.out_arb_intf.req.request;
   assign mout_arb_enable = s.out_arb_intf.req.enable;
   assign mout_arb_read = s.out_arb_intf.req.read;
   assign mout_arb_write = s.out_arb_intf.req.write;
   assign mout_arb_burstcount = s.out_arb_intf.req.burstcount;
   assign mout_arb_address = s.out_arb_intf.req.address;
   assign mout_arb_writedata = s.out_arb_intf.req.writedata;
   assign mout_arb_byteenable = s.out_arb_intf.req.byteenable;
   assign mout_arb_id = s.out_arb_intf.req.id;
   assign s.out_arb_intf.stall = mout_arb_stall;
   assign s.in_arb_intf.req = dp[0].out_intf.req;
   assign dp[0].out_intf.stall = s.in_arb_intf.stall;
   assign dp[0].in_intf.req = a[2].mout_intf.req;
   assign a[2].mout_intf.stall = dp[0].in_intf.stall;
   assign a[2].m0_intf.req = a[0].mout_intf.req;
   assign a[0].mout_intf.stall = a[2].m0_intf.stall;
   assign a[2].m1_intf.req = a[1].mout_intf.req;
   assign a[1].mout_intf.stall = a[2].m1_intf.stall;
   assign a[0].m0_intf.req = m[0].arb_intf.req;
   assign m[0].arb_intf.stall = a[0].m0_intf.stall;
   assign a[0].m1_intf.req = m[1].arb_intf.req;
   assign m[1].arb_intf.stall = a[0].m1_intf.stall;
   assign a[1].m0_intf.req = m[2].arb_intf.req;
   assign m[2].arb_intf.stall = a[1].m0_intf.stall;
   assign a[1].m1_intf.req = m[3].arb_intf.req;
   assign m[3].arb_intf.stall = a[1].m1_intf.stall;
endmodule

