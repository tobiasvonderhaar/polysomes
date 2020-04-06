def find_poly_peaks(trace,run_length=4, exclude_beginning=0.1):
    #read out some basic trace parameters
    trace_x_max = max(trace['x'])
    #ensure that x-values are ordered (this may not be the case with plots produced. via webplotdigitizer)
    trace = trace.sort_values('x')
    points = trace['y']
    points_prime = [points[y]-points[y-1] for y in range(1,len(points))]

    #determine positive runs
    pos_run_ends = []
    ploc = 0
    run_length=4
    while ploc < len(points_prime)-1:
        if points_prime[ploc] > 0:
            pos_start=ploc
            ploc+=1
            while points_prime[ploc] > 0 and ploc < len(points_prime)-1:
                ploc+=1
            pos_end = ploc
            if pos_end-pos_start > run_length:
                pos_run_ends.append(pos_end)
        else:
            ploc+=1

    #determine negative runs
    neg_run_starts = []
    nloc = 0
    while nloc < len(points_prime)-1:
        if points_prime[nloc] < 0:
            neg_start=nloc
            nloc+=1
            while points_prime[nloc] < 0 and nloc < len(points_prime)-1:
                nloc+=1
            neg_end = nloc
            if neg_end-neg_start > run_length:
                neg_run_starts.append(neg_start)
        else:
            nloc+=1

    #determine peak positions as the mean in x between a pos run end and the closest following neg run start
    global peakxs
    peakxs = []
    neg_start_ys = trace.iloc[neg_run_starts]['x'].values
    for entry in range(len(pos_run_ends)):
        pos_end_x = trace.iloc[pos_run_ends[entry]]['x']
        match_neg_start_x = min(neg_start_ys[neg_start_ys >= pos_end_x])
        #ensure that the pos run end and neg run start are not further than a maximum distance apart:
        if (match_neg_start_x - pos_end_x) < (trace_x_max / 30):
            #exclude cases where there is a second pos run end prior to the following matching neg run start
            #either include because this is the last value in pos_run_ends
            if entry == len(pos_run_ends)-1:
                print('ps_run_end length reached')
                peakxs.append((pos_end_x + match_neg_start_x)/2)
            #or include if there is no additional pos run between the current pos run and the next neg run
            elif match_neg_start_x < trace.iloc[pos_run_ends[entry + 1]]['x']:
                #exclude peaks near the beginning of the gradient, which are usually dirt
                if pos_end_x > trace_x_max * exclude_beginning:
                    peakxs.append((pos_end_x + match_neg_start_x)/2)
    return peakxs





def plot_trace(trace,peakxlocs=[],add_x=0,predxlocs=[]):
    import matplotlib.pyplot as plt
    fig,ax = plt.subplots()
    ax.scatter(trace['x'],trace['y'],c='black',s=2)
    if len(peakxlocs)>0:
        for loc in peakxlocs:
            ax.plot([loc,loc],[min(trace['y']),max(trace['y'])],c='green')
    if add_x > 0:
        ax.plot([add_x,add_x],[min(trace['y']),max(trace['y'])],c='orange')
    if len(predxlocs)>0:
        for loc in predxlocs:
            ax.plot([loc,loc],[min(trace['y']),max(trace['y'])],c='red')
    ax.set_ylabel(r'$oD_{254}$')
    plt.show()





def fit_peaks(peakxs,grad_end,mode='mammal'):
    from scipy.optimize import curve_fit
    import numpy as np
    #construct a dictionary of relative LSU/SSU MW data for different organisms
    LSU_SSU_data = {'mammal':0.37,'yeast':0.34}
    if type(mode)==str:
        LSU_SSU = LSU_SSU_data[mode]
    else:
        LSU_SSU = mode
    #construct a list with the molecular weights of consecutive polysomes (starting with 40S, 60S, monosome,disome,...)
    polyno = [LSU_SSU,1-LSU_SSU] + list(range(1,len(peakxs)-1))
    #define a logarithmic function to fit the peak x values to
    def logfunc(x,a,b,c):
        return a*np.log(b * x + c)
    #fit the parameters
    params, params_covariance = curve_fit(logfunc, polyno, peakxs)
    #calculate the fitted position of peaks following the observed peaks
    xfit=np.arange(max(polyno)+1,50)
    peakfit = logfunc(xfit,params[0],params[1],params[2])
    peakfit = peakfit[peakfit<grad_end]
    xfun=np.linspace(0.33,len(peakxs)+len(peakfit),150)
    peakfun = logfunc(xfun,params[0],params[1],params[2])
    return [peakfit, [xfun,peakfun]]

def fp2poly(dataset, from_file = True, use_ref_RNA = False, RNA_content = 60000, ribo_content = 200000, frac_act = 0.85, frac_split = 0.30, poly_limit=30, remove_spurious_RNAs = True,):
    """Calculates peak volumes of a polysome profile corresponding to an input footprinting dataset."""
    
    import pandas as pd
    import numpy as np
    
    genes = pd.read_csv('Data/sacCer3 genes.csv')
    
    if from_file:
        ######read in dataset and validate that it has columns with required names
        if dataset[-4:] != '.csv':
            filename = dataset + '.csv'
        else:
            filename = dataset
        dats = pd.read_csv("Data/" + filename)
    else:
        dats = dataset
        
    if ('ORF' and 'Ribo_Prints') not in dats.columns:
        print('One or more required columns are missing.')
        return
    
    if use_ref_RNA:
        if 'RNA_prints' in dats.columns:
            dats = dats.drop('RNA_Prints',axis=1)
    
    if 'RNA_Prints' not in dats.columns:
        RNA_ref = pd.read_csv('Data/RNA_reference.csv')
        dats = dats.merge(RNA_ref,how='inner',on='ORF')
        print('Using reference mRNA data.')
    
    #remove datasets where either the RNA prints or Ribo prints are 0
    dats = dats.loc[dats['RNA_Prints'] > 0]
    dats = dats.loc[dats['Ribo_Prints'] > 0]
    
    ######convert reads to RPKM
    
    #combine input dataset with gene length information 
    dats = dats.merge(genes[['name','length']],how='inner',left_on='ORF',right_on='name')[['ORF','RNA_Prints','Ribo_Prints','length']]
    #calculate RPKM
    dats['RNA_RPKM'] = dats['RNA_Prints']/(dats['length']/1000)    
    
    #####sort genes into polysome peaks according to RPKM info
    
    #determine conversion factor from Ribo_Prints to no of Ribosomes
    RiboPrints2Ribos = (ribo_content * frac_act) / sum(dats['Ribo_Prints'])
    dats['Ribos_bound'] = dats['Ribo_Prints']*RiboPrints2Ribos
    #determine conversion factor from RNA_RPKM to no of RNAs
    RNARPKM2RNAs = RNA_content / sum(dats['RNA_RPKM'])
    dats['RNAs_per_cell'] = dats['RNA_RPKM']*RNARPKM2RNAs
    #calculate the ribosome load per RNA (RperR)
    dats['RperR'] = dats['Ribos_bound'] / dats['RNAs_per_cell']
    dats=dats.dropna()
    #remove rows where the number of ribosomes per RNA is > poly_limit
    dats = dats.loc[dats['RperR'] <= poly_limit]
    #remove spurious RNAs (< than 0.05 RNAs per cell)
    if remove_spurious_RNAs:
        dats = dats.loc[dats['RNAs_per_cell'] > 0.05]
    
    ######assign RNAs into polysome peaks
    
    #make an array to hold the relative weights for each polysome class
    poly_array = np.zeros(poly_limit+2)
    #assign idle ribosomes to the first three peaks, based on the fraction of active ribosomes and the fraction of split inactive ribosomes
    idle_ribos = (1 - frac_act) * 200000
    poly_array[0] += idle_ribos * frac_split * 0.34
    poly_array[1] += idle_ribos * frac_split * 0.66
    poly_array[2] += idle_ribos * (1-frac_split)
    #go through each row of dats and add ribosomes to the appropriate peak
    for row in range(dats.shape[0]):
        this_RperR = dats.iloc[row]['RperR']
        these_Ribos_bound = dats.iloc[row]['Ribos_bound']
        #if the number of ribos per RNA is an exact number, assign the ribos to the corrresponding peak
        floor_val = int(this_RperR)
        if  float(floor_val) == this_RperR:
            poly_array[floor_val + 1] += these_Ribos_bound
        #if the number of ribos is between two integers, split the ribos proportionally between the two adjacent peaks
        #for example, for 5.6 Ribos per RNA 60% of ribosomes go to the 6-some peak, 40% to the 5-some peak
        else:
            ceil_weight = (this_RperR-floor_val)
            floor_weight = 1 - ceil_weight
            if floor_val != 0:
                poly_array[floor_val + 1] += these_Ribos_bound * floor_weight
            poly_array[floor_val + 2] += these_Ribos_bound * ceil_weight
    #normalise values to total signal
    poly_array = poly_array/sum(poly_array)
    return poly_array

def plot_poly(peak_vols):
    """Plots a polysome profile from a list of peak volumes computed by fp2poly"""
    
    import numpy as np
    
    #define the function for returning a normal distribution centred around mu with variance (=peak width) sigma
    def normpdf(x,mu,sigma):
        return 1/(np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(-1 * ((x - mu) ** 2 / (2 * sigma ** 2)))
    
    #define the relative molecular weights corresponding to the different peaks (40S, 60S, 1-,2-,3-,... some)
    peak_ids = [0.34,0.66] + list(range(1,len(peak_vols)-1))
    #calculate a series of peak locations based on a typical polysome profile
    peak_locs = [0.15*np.log(19 * x -0.6) for x in peak_ids]
    #calculate a series of peak widths based on a typical polysome profile
    peak_widths = [-0.0004 * x + 0.015 for x in peak_ids]
    #definine the oD drift and the initial peak based on a typical polysome profile
    x = np.linspace(0,1,num=400)
    drift = -0.07 + 0.19*x + 0.05
    debris = np.exp((-x+0.085)*20)*1.2
    #construct the plot of the polysome profile
    sum_trace = drift + debris
    for peak_no in range(len(peak_locs)):
        this_peak = normpdf(x,peak_locs[peak_no],peak_widths[peak_no])*peak_vols[peak_no]
        sum_trace += this_peak
    return x,sum_trace
    
    
    
    
def compare_profiles(dats1, dats2, dats1_columns = ['ORF','RNA_Prints','Ribo_Prints'],
                     dats2_columns = ['ORF','RNA_Prints','Ribo_Prints'],colors =['royalblue','gold'],conditions = ['Cond. 1','Cond. 2'],return_df=False):
    
    """Computes and displays the predominant movements of transcripts between polysome peaks for two conditions."""
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    def counts_to_RPKM(dats):
        #read in the gene length info dataset
        genes = pd.read_csv('Data/sacCer3 genes.csv')
        #combine input dataset with gene length information 
        dats = dats.merge(genes[['name','length']],how='inner',left_on='ORF',right_on='name')[['ORF','RNA_Prints','Ribo_Prints','length']]
        #calculate RPKM
        dats['RNA_RPKM'] = dats['RNA_Prints']/(dats['length']/1000) 
        dats = dats.drop(['RNA_Prints','length'],axis=1)
        return dats
    
    def calc_RperR(Rdats, RNA_content = 60000, ribo_content = 200000, frac_act = 0.85,poly_limit = 30):
        #determine conversion factor from Ribo_Prints to no of Ribosomes
        RiboPrints2Ribos = (ribo_content * frac_act) / sum(Rdats['Ribo_Prints'])
        Rdats['Ribos_bound'] = Rdats['Ribo_Prints']*RiboPrints2Ribos
        #determine conversion factor from RNA_RPKM to no of RNAs
        RNARPKM2RNAs = RNA_content / sum(Rdats['RNA_RPKM'])
        Rdats['RNAs_per_cell'] = Rdats['RNA_RPKM']*RNARPKM2RNAs
        #calculate the ribosome load per RNA (RperR)
        Rdats['RperR'] = np.round(Rdats['Ribos_bound'] / Rdats['RNAs_per_cell'])
        Rdats=Rdats.dropna()
        #remove rows where the number of ribosomes per RNA is > poly_limit
        Rdats = Rdats.loc[Rdats['RperR'] <= poly_limit]
        return Rdats
    
    #prepare datasets and compute RperR values via RNA RPKM values
    dats1 = dats1[dats1_columns]
    dats1.columns = ['ORF','RNA_Prints','Ribo_Prints']
    dats1 = counts_to_RPKM(dats1) 
    dats1 = calc_RperR(dats1)
    dats2 = dats2[dats2_columns]
    dats2.columns = ['ORF','RNA_Prints','Ribo_Prints']
    dats2 = counts_to_RPKM(dats2)
    dats2 = calc_RperR(dats2)
    
    #compute column movements for individual transcripts
    fromto = dats1[['ORF','RperR']].merge(dats2[['ORF','RperR']],how = 'inner',on = 'ORF')
    fromto.columns = ['ORF','from','to']
    
    #calculate main destinations for each origin peak
    from_unique = fromto['from'].unique()
    from_unique.sort()
    from_vec,to_vec,dir_vec,alpha_vec=[],[],[],[]
    for from_idx in range(len(from_unique)):
        this_from  = fromto.loc[fromto['from'] == from_unique[from_idx]]
        this_to_unique = this_from['to'].unique()
        for to_idx in range(len(this_to_unique)):
            from_vec.append(from_unique[from_idx])
            to_vec.append(this_to_unique[to_idx])
            if from_unique[from_idx] < this_to_unique[to_idx]:
                dir_vec.append('up')
            elif from_unique[from_idx] > this_to_unique[to_idx]:
                dir_vec.append('down')
            else:
                dir_vec.append('nc')
            alpha_vec.append(np.log(this_from.loc[this_from['to']==this_to_unique[to_idx]].shape[0]))
    alpha_vec_scaled = [(alpha)/(max(alpha_vec)*0.7) for alpha in alpha_vec]
            
    df = pd.DataFrame({'From':from_vec,'To':to_vec,'Direction':dir_vec,'Alpha':alpha_vec_scaled})
    up_df = df.loc[df['Direction']=='up']
    down_df = df.loc[df['Direction']=='down']
     
    #prepare main figure
    fig = plt.figure(constrained_layout=True,figsize=(5,3))
    gs = fig.add_gridspec(2, 6)
    ax1 = fig.add_subplot(gs[0, :-2])
    ax2 = fig.add_subplot(gs[1, :-2])
    ax3 = fig.add_subplot(gs[0:, -2:])
    ax3.axis('off')
    
    for idx in range(up_df.shape[0]):
        ax1.plot([up_df.iloc[idx]['From'],up_df.iloc[idx]['To']],[3,1],linewidth=4,color=colors[0],alpha=up_df.iloc[idx]['Alpha'])
    ax1.set_ylim(0.9,3.1)
    ax1.set_yticks([3,1])
    ax1.set_yticklabels([conditions[0],conditions[1]])
    ax1.set_xlim((0,30))
    
                                                                                         
    for idx in range(down_df.shape[0]):
        ax2.plot([down_df.iloc[idx]['From'],down_df.iloc[idx]['To']],[3,1],linewidth=4,color=colors[1],alpha=down_df.iloc[idx]['Alpha'])
    ax2.set_ylim(0.9,3.1)
    ax2.set_yticks([3,1])
    ax2.set_yticklabels([conditions[0],conditions[1]])
    ax2.set_xlabel('Polysome number')
    ax2.set_xlim((0,30))
    
    #prepare legend for line width
    ax3.set_ylim((0,10))
    ax3.set_xlim((0,10))
    ax3.text(8.5,5,'Number of mRNAs',verticalalignment='center',rotation = 90)
    ref_alphas = [0.1,0.4,0.7]
    ref_values = [np.exp(max(alpha_vec)/10),np.exp(max(alpha_vec)/2),np.exp(max(alpha_vec))]
        
    y_pos = [3.5,5,6.5]
    for idx in [0,1,2]:
        ax3.plot([0.5,3.5],[y_pos[idx],y_pos[idx]],linewidth=4,c='black',alpha = ref_alphas[idx])
        ax3.text(4.5,y_pos[idx],str(int(ref_values[idx])),verticalalignment='center')
    
    if return_df:
        return fig,fig.axes,df
    else:
        return fig,fig.axes