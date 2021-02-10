import csv
import pandas as pd
import numpy as np
#
def find_target(regs, ext_cut_off):

def cluster_reg_element(regs, loops, ext_cut_off):
  """
  Function: Cluster regulatory elements (enhancers, promoters) by 
            the structural constraints of chromatin loops
  Inputs:
    regs: a dataframe, each row represents a regulatory element, 
          3 columns "chr", "start", "end" represent the coordinate
    loops: a dataframe, each row represents a loop, 
           3 columns "chr", "start", "end" represent the coordinate
    ext_cut_off: a threshold value (number of base pairs)
  Output:
    a list of lists, each list contains row ids of the input 
    dataframe regs to represent a cluster of regulatory elements 
    of these row ids
  """
  regs["origin_index"] = regs.index
  clusters = []
  for chr in regs.chr.unique():
    #print(chr)
    chr_spec_regs = regs[(regs['chr'] == chr)]
    chr_spec_regs.reset_index(inplace=True)
    chr_spec_loops = loops[(loops['chr'] == chr)]
    # 100kb only appears once here
    chr_spec_shortest_ext = find_shortest_ext_len(chr_spec_regs, chr_spec_loops, 100000)
    #print (chr_spec_shortest_ext)
  return cluster
def find_shortest_ext_len(regs, loops, max_ext_len):
  """
  Function: Find the shortest extension length so that a pair of 
            regulatory elements can interact  
  Inputs:
    regs: a dataframe, the column "start" represent the position 
          of a regulatory element represented by a row
    loops: a dataframe, two columns "start" and "end" represent 
           the coordinate of a loop represented by a row
    max_extension: a threshold value (number of base pairs) 
  Output:
    a dataframe, 3 columns "row_1", "row_2", "shortest_ext_len"
  
  """
  reg_num = len(regs.index)
  min_loop_in = np.full((reg_num, reg_num), 1e7)
  max_loop_out = np.zeros((reg_num, reg_num))
  # add
  for loop_index,loop in loops.iterrows():
    for i, reg_i in regs.iterrows():
      for j, reg_j in regs.iterrows():
        if (i > j):
          min_loop_in[i][j] = min(min_loop_in[i][j], 
                                  find_loop_in_extension(loop['start'], loop['end'], 
                                                         reg_i['start'], reg_j['start']))
  # remove
  for loop_index,loop in loops.iterrows():
    tmp = []
    for i,reg_i in regs.iterrows():
      for j,reg_j in regs.iterrows():
        if (i > j):
          max_loop_out[i][j] = max(max_loop_out[i][j],
                                   find_loop_out_extension(loop['start'], loop['end'], 
                                                           reg_i['start'], reg_j['start']))
  # combine
  tmp = []
  for j in range(reg_num):
    for i in range(j + 1, reg_num):
      shortest_ext_len = max(min_loop_in[i][j], max_loop_out[i][j])
      if (shortest_ext_len <= max_ext_len):
        tmp.append([i, j, shortest_ext_len])
  return pd.DataFrame.from_records(tmp, columns=["row_1", "row_2", "shortest_ext_len"])

# 
def find_loop_in_extension(a1, a2, e1, e2):
  if (e1 >= a2):
    return e2 - a2
  elif (e1 >= a1):
    if (e2 >= a2):
      return e2 - a2
    else:
      return 0
  else:
    if (e2 <= a1):
      return a1 - e1
    elif (e2 <= a2):
      return a1 - e1
    else:
      return max(a1 - e1, e2 - a2)
#
def find_loop_out_extension(a1, a2, e1, e2):
  if (e1 >= a1 and e1 <= a2 and e2 > a2):
    return e2 - a2
  elif (e1 < a1 and e2 >= a1 and e2 <= a2):
    return a1 - e1
  else:
    return 0

## Tests
def test(test_num):
  loop_df = pd.read_csv("../exp/2021_01_05_dot_call/dot_mic_mESC_5kb.postproc", sep='\t')
  loop_df["chr"] = loop_df["chrom1"]
  loop_df["start"] = round((loop_df["start1"] + loop_df["end1"])/2)
  loop_df["end"] = round((loop_df["start2"] + loop_df["end2"])/2)
  loops = loop_df[["chr", "start", "end"]]
  print(loops)
  # simple test
  if (test_num >= 1):
    reg_short_list_df = pd.read_csv("../data/2021_01_14_E_P_mESC/Jennifer_list.csv", sep = '\t')
    print(cluster_reg_element(reg_short_list_df, loops, 5000))

  # long test
  if (test_num >= 2):
    enhancer_df = pd.read_csv("../../data/2021_01_18_Jennifer_enhancers/Jennifer_Supplemental_Table_S8.csv", sep='\t')
    enhancer_df.rename(columns = {"Chromosome": "chr"}, inplace=True)
    enhancer_df["id"] = "enh_" + enhancer_df.index.astype(str)
    enhancer_df["type"] = "enhancer"
    promoter_df = pd.read_csv("../exp/2021_01_18_biomart/results/principal_mm10_ensemble_genes.csv", sep=' ')
    promoter_df.rename(columns = {"chromosome_name": "chr", "external_gene_name": "id", "transcription_start_site": "start"}, inplace=True)
    promoter_df['chr'] = 'chr' + promoter_df['chr'].astype(str)
    promoter_df['end'] = promoters['start'].astype(int) + 100
    promoter_df['type'] = 'promoter'
    reg_long_list_df = pd.concat([enhancer_df[['id', 'chr', 'start', 'end', 'type']], promoter_df[['id', 'chr', 'start', 'end', 'type']]])
    #print(reg_element_long_list)
    #reg_element_long_list = reg_element_long_list[(reg_element_long_list['chr'] == 'chr10')]
    reg_long_list_df.reset_index(drop=True, inplace=True)
    print(cluster_reg_element(reg_short_list_df, loops, 5000))
##
test(1)
