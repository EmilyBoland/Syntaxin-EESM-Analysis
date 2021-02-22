def calc_work(alldata, mem, run):
  """
  Calculates the work for one simulation run of one member. 
  inputs: 
          alldata = json with all data,
          mem = the member number, 
          run = the simulation run 
  output: work for all three distributions of one run of one member
  """

  work = 0 
  dimensions = ['052_210', '105_216', '196_228']     
  for dim in dimensions:
    Rs = alldata['converged run data'][mem][run][dim]['R']
    distances = []
    for i in range(0,len(Rs)-1):
      lowerR = Rs[i]
      upperR = Rs[i+1]
      distance = upperR-lowerR
      distances.append(distance)
        
    alpha = alldata['converged run data'][mem][run][dim]['alpha']
    target = alldata['converged run data'][mem][run][dim]['target']
    force = alpha/target 
    for distance in distances:
      work += abs(distance*force)
        
  return work    




def threshold(threshold, work_values, work_values_by_run):
  """ sets all work values less than the threshold value equal to the new threshold value
    inputs: 
      threshold = value to use as threshold (int)
      work_values = dict containing target_set: list of work values pairs
      work_values_by_run = nested dict containing target_set: list of work values pairs for each run
    output: nested dict containing updated work_values and work_values_by_run
  """
  for target_set in work_values: 
    change_ws = []
    for w in work_values[target_set]:
      if w < threshold: 
        change_ws.append(w)
    for change_w in change_ws: 
      work_values[target_set].remove(change_w)
      work_values[target_set].append(threshold)
    for run in work_values_by_run: 
      try:
        for w in work_values_by_run[run][target_set]:
          if w < threshold: 
            work_values_by_run[run][target_set].remove(w)
            work_values_by_run[run][target_set].append(threshold)
      except:
        pass
  return {'work_values': work_values, 'work_values_by_run': work_values_by_run}




def average_jarz_work(w_vals):
  """
  implements the Jarzynski Equation for one target 
  inputs: 
    w_vals = list of work values for a single set of targets 
  output: exponential work value (e^(-beta*work)) averaged across all values in the input list
  """

  import math
  
  # beta = 1/(R*T), where R is the ideal gas constant in kJ/(mol*K) and T is standard human body temperature: 310 K
  R = 0.0083144598
  T = 310 
  beta = 1/(R*T)
  
  n = 0
  running_w_total = 0
  
  if len(w_vals) == 0: 
    avg_jarz_w = 0
  
  else:
    for w in w_vals: 
      if w == 'INF':
        n += 1
      else:
        w = float(w)
        running_w_total += math.exp(-beta*w)
        n += 1
    avg_jarz_w = running_w_total/n

  return avg_jarz_w




def normalize_probs(jarz_work_dict):
  """
  normalizes probabilities of each configuration (for each set of targets)
  input: 
    [if dictionary]   jarz_work_dict = dictionary containing average exponential work terms given by the Jarzynski Equation value for each target (see above)
    [if list]         list of liklihoods that correlate to probabilities when normalized to 1
  output: 
    [if dictionary]   normalized_probs = dictionary containing probability of a given target: {'targets': probability}
    [if list]         list that sums to one so that all items within the list are normalized
  """
  if type(jarz_work_dict) is list: 
    normalized_probs = []
    for item in jarz_work_dict: 
      normalized_probs.append(item/sum(jarz_work_dict))
  else:
    Psum = sum(list(jarz_work_dict.values()))
    normalized_probs = {}
    for key in jarz_work_dict:
      normalized_probs[key] = jarz_work_dict[key]/Psum
  return normalized_probs




def twoD_sum(data): 
  """ sums over all 3 dimensions in a 3D array to make 3x 2D arrays
    input: 
      data = a 3D array (same dimensions in all 3 directions)
    returns: 
      a dict containing 3x 2D arrays obtained 
  """
  
  import numpy as np
  
  sum_dict = {}
  sum_dict['052_210__105_216'] = np.sum(data, axis = 2)
  sum_dict['052_210__196_228'] = np.sum(data, axis = 1)
  sum_dict['105_216__196_228'] = np.sum(data, axis = 0) 
  
  return sum_dict




def make_map(plot_dict, colorchoice = "hot"):
  """
  plots the data given in a 2-dimensional array
  inputs: 
    array = 2-D array containing data to be plotted
    xlabel = label for x axis
    ylabel = label for y axis
    norm = Normalization: defaults to none and assumes linear data, can be maniputated for log scale
  output: the boolean True, in running, it plots the data
  """

  import matplotlib.pyplot as plt
  import matplotlib.gridspec as gridspec
  import numpy as np

  fig = plt.figure(figsize = (20,20))
  spec = gridspec.GridSpec(ncols=3, nrows=1, figure = fig, wspace = 0.5)
  i = 1
  for key in plot_dict:
    xlabel = key[9:]
    ylabel = key[0:7]
    if colorchoice == 'hot':
      fig.add_subplot(3,3,i)
      plt.imshow(plot_dict[key], 
        interpolation=None, 
        origin='lower', cmap= colorchoice)
      plt.xlabel(xlabel)
      plt.ylabel(ylabel)
      plt.colorbar()
      plt.title("{} v. {}".format(ylabel, xlabel))
      plt.xticks(np.arange(16), ('','1.0', '','2.0', '','3.0','', '4.0','', '5.0','','6.0','','7.0','','8.0'))
      plt.yticks(np.arange(16), ('','1.0', '','2.0', '','3.0','', '4.0','', '5.0','','6.0','','7.0','','8.0'))
    else:
      bound = max(max(abs(x)) for x in plot_dict[key])
      xlabel = key[9:]
      ylabel = key[0:7]
      fig.add_subplot(3,3,i)
      plt.imshow(plot_dict[key], 
        interpolation=None, 
        origin='lower', cmap= colorchoice, 
        vmin = -bound, vmax = bound)
      plt.xlabel(xlabel)
      plt.ylabel(ylabel)
      plt.colorbar()
      plt.title("{} v. {}".format(ylabel, xlabel))
      plt.xticks(np.arange(16), ('','1.0', '','2.0', '','3.0','', '4.0','', '5.0','','6.0','','7.0','','8.0'))
      plt.yticks(np.arange(16), ('','1.0', '','2.0', '','3.0','', '4.0','', '5.0','','6.0','','7.0','','8.0'))
    i += 1 
  return True




def oneD_sum(data): 
  """ sums over all 3 dimensions in a 3D array to make 3x 2D arrays
    input: 
      data = a 3D array (same dimensions in all 3 directions)
    returns: 
      a dict containing 3x 2D arrays obtained 
  """
  
  import numpy as np
  
  sum_dict = {}
  sum_dict['052_210'] = np.sum(data, axis = (1,2))
  sum_dict['105_216'] = np.sum(data, axis = (0,2))
  sum_dict['196_228'] = np.sum(data, axis = (0,1)) 
  
  return sum_dict




def JScalculator(pi, pf):
  """
  calculates JS divergence between two given distributions
  inputs: 
    pf = a list of probabilities of length n corresponding to n bins
    pi = a second list of probabilities of length n corresponding to n bins
  outputs: 
    a float indicating degree of similarity between pf and pi
  """

  from scipy.stats import entropy
  
  ms = []
  for i in range(len(pi)):
    ms.append(0.5 * pi[i] + 0.5 * pf[i])
  jsdiv = 0.5 * (entropy(pi, qk = ms) + entropy(pf, qk = ms))
  
  return jsdiv




def plot_js_along_runs(by_run_dict, final_dict, final_dict_1D, dat_type):
  """
  creates a plot of JS divergence between final distribution v. simulation data along runs (cumulative)
  inputs: 
    by_run_dict = the data by run (cumulative) to calculate JS divergence with
    final_dict = the final data to calculate JS divergence with 
    final_dict_1D = a 3 key:value pair dictionary containing 3 lists with the final data summed along all three axes
    dat_type = str describing source of data in final_dict. Either "Experimental" or "Simulation"
  output: 
    a plot containing 4 curves: one for each distribution and an overall
  """

  import matplotlib.pyplot as plt
  import numpy as np

  by_run_1D = {'052_210':{}, '105_216': {}, '196_228':{}}
  for run in by_run_dict:
    by_run_1D['052_210'][run] = np.sum(by_run_dict[run], axis =(1,2))
    by_run_1D['105_216'][run] = np.sum(by_run_dict[run], axis =(0,2))
    by_run_1D['196_228'][run] = np.sum(by_run_dict[run], axis =(0,1))

  total_js_divs = []
  by_run_lst_dict = {}
  for run in by_run_dict:
    final_lst = []
    by_run_lst = []
    by_run_data = by_run_dict[run]
    for d1 in range(17):
      for d2 in range(17):
        for d3 in range(17):
          final_lst.append(final_dict[d1,d2,d3])
          by_run_lst.append(by_run_data[d1,d2,d3])
    total_js_divs.append(JScalculator(final_lst, by_run_lst))

  jsdivs_by_run = {}
  for distr in by_run_1D:
    jsdivs_by_run[distr] = []
    for run in by_run_1D[distr]:
      jsdivs_by_run[distr].append(JScalculator(by_run_1D[distr][run], final_dict_1D[distr]))

  runs = [1,2,3,4,5,6,7,8,9,10]
  plt.plot(runs,jsdivs_by_run['052_210'],':', c = "purple", lw = 3)
  plt.plot(runs,jsdivs_by_run['105_216'],'--c', lw = 3)
  plt.plot(runs,jsdivs_by_run['196_228'],'-.', c = 'orange', lw = 3)
  plt.plot(runs, total_js_divs, '-r', lw = 3)
  plt.xlabel("Number of Runs (Cumulative)")
  plt.ylabel("JS-Divergence")
  plt.title("JS-Divergence between final {} Joint Distribution \n and Simulation Distribution along Number of Runs (Cumulative)".format(dat_type))
  plt.legend(['052_210', '105_216', '196_228', 'total'])

