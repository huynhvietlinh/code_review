import csv
import pandas as pd
import numpy as np
#
def find_target(reg_ids, interactions, ext_cut_off):
  """
  Function:
  Inputs:
    reg_ids:
    interactions:
    ext_cut_off: a threshold value (number of base pairs)
  Output:
    a list of ids
  """
  targets = []
  return targets
#
def find_interactions(regs, loops, max_ext_len = 1e6):
  """
  Function: find possible interactions between regulatory elements (enhancers, promoters) by 
            the structural constraints of chromatin loops
  Inputs:
    regs: a dataframe, each row represents a regulatory element, 
          3 columns "chr", "start", "end" represent the coordinate
    loops: a dataframe, each row represents a loop, 
           3 columns "chr", "start", "end" represent the coordinate
    max_ext_len: a threshold value (number of base pairs) 
  Output:
    a dataframe, each row represents a possible interaction,
    2 columns "row_1", "row_2" are row indexes of two regulatory 
    elements in the input dataframe regs and column "min_ext_len" is the 
    minimum extension length so that these two elements can interact
  """
  
  regs["origin_index"] = regs.index
  interaction_list = []
  for chr in regs.chr.unique():
    #print(chr)
    chr_loops = loops[(loops["chr"] == chr)]
    chr_regs = regs[(regs["chr"] == chr)]
    # sort "start" increasingly
    chr_regs = chr_regs.sort_values("start")
    chr_regs.reset_index(drop=True, inplace=True)
    chr_reg_num = len(chr_regs.index)
    min_loop_in = np.full((chr_reg_num, chr_reg_num), 1e7)
    max_loop_out = np.zeros((chr_reg_num, chr_reg_num))
    # elements in a loop can interact
    for loop_index,loop in chr_loops.iterrows():
      for i, reg_i in chr_regs.iterrows():
        for j, reg_j in chr_regs.iterrows():
          if (i < j):
            min_loop_in[i][j] = min(min_loop_in[i][j], 
                                    find_loop_in_extension(loop["start"], loop["end"], 
                                                           reg_i["start"], reg_j["start"]))
    # elements separated by a loop can not interact
    for loop_index,loop in chr_loops.iterrows():
      for i,reg_i in chr_regs.iterrows():
        for j,reg_j in chr_regs.iterrows():
          if (i < j):
            max_loop_out[i][j] = max(max_loop_out[i][j],
                                     find_loop_out_extension(loop["start"], loop["end"], 
                                                             reg_i["start"], reg_j["start"]))
    # combine
    for i in range(chr_reg_num):
      for j in range(i + 1, chr_reg_num):
        shortest_ext_len = max(min_loop_in[i][j], max_loop_out[i][j])
        if (shortest_ext_len <= max_ext_len):
          interaction_list.append([chr_regs.iloc[i]["origin_index"], chr_regs.iloc[j]["origin_index"], 
                                   shortest_ext_len])
  #
  return pd.DataFrame.from_records(interaction_list, columns=["row_1", "row_2", "shortest_ext_len"])
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
def find_loop_out_extension(a1, a2, e1, e2):
  if (e1 >= a1 and e1 <= a2 and e2 > a2):
    return e2 - a2
  elif (e1 < a1 and e2 >= a1 and e2 <= a2):
    return a1 - e1
  else:
    return 0

 
def find_loop_in_extension(LS, LE, e1, e2):
  """
  LS, LE: loop start and loop end
  e1, e2: start position of two elements, e1 <= e2
  """
  if (e1 >= LE):
    return e2 - LE #---LS---LE---e1---e2----
  elif (e1 >= LS):
    if (e2 >= LE):
      return e2 - LE #---LS---e1---LE---e2----
    else:
      return 0 #---LS----e1--e2----LE----
  else:
    if (e2 <= LS):
      return LS - e1 #---e1---e2---LS----LE----
    elif (e2 <= LE):
      return LS - e1 #---e1---LS----e2----LE----
    else:
      return max(LS - e1, e2 - LE) #---e1---LS---LE----e2-----
#
def find_loop_out_extension(LS, LE, e1, e2):
  """
  LS, LE: loop start and loop end
  e1, e2: start position of two elements, e1 <= e2
  """
  if (e1 >= LS and e1 <= LE and e2 > LE):
    return e2 - LE #---LS---e1----LE----e2----
  elif (e1 < LS and e2 >= LS and e2 <= LE):
    return LS - e1 #---e1----LS----e2----LE---
  else:
    return 0

## Tests
def test(test_num):
  loop_df = pd.read_csv("../exp/2021_01_05_dot_call/dot_mic_mESC_5kb.postproc", sep='\t')
  loop_df["chr"] = loop_df["chrom1"]
  loop_df["start"] = round((loop_df["start1"] + loop_df["end1"])/2)
  loop_df["end"] = round((loop_df["start2"] + loop_df["end2"])/2)
  loops = loop_df[["chr", "start", "end"]]
  #print(loops)
  # simple test
  if (test_num >= 1):
    reg_short_list_df = pd.read_csv("../data/2021_01_14_E_P_mESC/Jennifer_list.csv", sep = '\t')
    print(find_interactions(reg_short_list_df, loops))

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
    print(find_interactions(reg_short_list_df, loops))
##
test(1)
