import pandas as pd
import numpy as np
from scipy import stats
import streamlit as st
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool, Range1d
##########################################################################################################################################
class Experiment:
    
    def __init__(self):
        self.n_conditions = None
        self.conditions_list = []
        self.number_of_samples_list = []
        self.individual_samples_list = []
        self.full_samples_list = []
        
    def get_input(self):
        
        st.sidebar.subheader("Define Experiment")
        
        # Number of experimental groups or conditions
        def get_n_conditions():
            return st.sidebar.number_input('Enter the number of conditions',min_value = 1, max_value= 20, value = 1, step = 1)
        
        self.n_conditions = get_n_conditions()

        # List of conditions
        def get_conditions_list():
            return [st.sidebar.text_input('Create a label for condition #'+str(i)+' (e.g. WT or KO)') for i in range(self.n_conditions)]
        
        self.conditions_list = get_conditions_list()
        
        # print an error message if a condition has no label
        if "" in self.conditions_list:
            raise Exception ("Condition's label must be at least one character long! ")

        # List of the number of samples corresponding to each condition
        def get_number_of_samples_list():
            return [st.sidebar.number_input('Enter the number of samples for condition #'+str(i), min_value = 1, max_value = 1000, value = 1, step = 1)\
                    for i in range(self.n_conditions)]
        
        self.number_of_samples_list = get_number_of_samples_list()
        
        # List of all samples 
        def create_full_samples_list():
            return ['s' + str(i+1) for i in range(sum(self.number_of_samples_list))]
        
        self.full_samples_list = create_full_samples_list()
        
        # List of the list of samples corresponding to each condition 
        def create_individual_samples_list():
            temporary_list = self.full_samples_list.copy() # Make a copy of full_sample_list to prevent affecting it 
            for condition, replicate in zip(self.conditions_list, self.number_of_samples_list):
                self.individual_samples_list.append(temporary_list[: replicate])
                del temporary_list[: replicate]
            return self.individual_samples_list
        
        self.individual_samples_list = create_individual_samples_list()
            
    def __repr__(self):
        return " conditions_list:{},\n number_of_samples_list:{},\n individual_samples_list:{},\n full_samples_list:{},\n"\
        .format(self.conditions_list, self.number_of_samples_list, self.individual_samples_list, self.full_samples_list)
    
    def __str__(self):
        return self.__repr__()
    
###########################################################################################################################
class CleanData:
    
    def __init__(self, dataframe, experiment_object):
        
        self.dataframe = dataframe
        
        self.experiment_object = experiment_object 
        
        # Check the validity of user inputs
        # Compare the total number of 'MainArea' columns in df with total number of replicates given by the user 
        if len([col for col in self.dataframe.columns if 'MainArea[s' in col]) != \
        len(self.experiment_object.full_samples_list):
            raise Exception ("Invalid total number of samples!")
        
    def clean_data(self):

        # Extract relevant columns for analysis 
        def extract_relevant_columns():
            return self.dataframe[['Rej','LipidMolec', 'Class', 'Calc Mass', 'BaseRt'] + 
                                  ['MainArea[' + sample + ']' for sample in self.experiment_object.full_samples_list]]
        
        clean_dataframe = extract_relevant_columns()
        
        # Apply filter
        def apply_filter():
            # removes the datapoint if 'Rej' = 1
            return clean_dataframe.loc[clean_dataframe['Rej'] == 0]
        
        clean_dataframe = apply_filter()
        
        # Impute missing values: replace 0 values by the smallest non-zero value in the column 
        def impute_missing_value(clean_dataframe):
            
            for sample in self.experiment_object.full_samples_list:

                non_zero_list = [ele for ele in clean_dataframe['MainArea[' + sample + ']'].values if ele > 0]

                impute_value = min(non_zero_list)

                clean_dataframe['MainArea[' + sample + ']'] = \
                clean_dataframe['MainArea[' + sample + ']'].apply(lambda x: impute_value if x<=0 else x)

            return clean_dataframe

        clean_dataframe = impute_missing_value(clean_dataframe)

        return clean_dataframe
    
###########################################################################################################################
class VolcanoPlot:
    
    def __init__(self, experiment_object):
        self.control = st.selectbox('pick the control condition', experiment_object.conditions_list, 0)
        self.experimental = st.selectbox('Pick the experimental condition', experiment_object.conditions_list, 1)
        self.p_value_threshold = st.number_input('Enter the significance level', min_value = 0.001, max_value= 0.1, value = 0.01, step = 0.001)
        self.q_value_threshold = -np.log10(self.p_value_threshold)
        
    def add_fold_change_and_p_value_columns(self, clean_dataframe, experiment_object):
        
        def fold_change_calculator(nums_1, nums_2): 
            return np.log2(np.mean(nums_1)/np.mean(nums_2))
        
        def p_value_calculator(nums_1, nums_2):
            t_value, p_value = stats.ttest_ind(nums_1,nums_2)
            return p_value
        
        # get index of control and experimental conditions in experiment_object.conditions_list
        def get_index():
            return experiment_object.conditions_list.index(self.control),\
                   experiment_object.conditions_list.index(self.experimental)
            
        control_idx, experimental_idx = get_index()
        
        # get list of control samples and experimental samples from experiment_object.individual_samples_list
        def get_individual_samples_list():
            return experiment_object.individual_samples_list[control_idx], \
                   experiment_object.individual_samples_list[experimental_idx]
        
        control_samples_list, experimental_samples_list = get_individual_samples_list()
                        
        # add fold change column to clean_dataframe using "fold_change_calculator"
        clean_dataframe['fc_' + self.experimental + '_' + self.control] = \
        clean_dataframe[['MainArea[' + sample + ']'  for sample in control_samples_list + experimental_samples_list]]\
        .apply(lambda x: fold_change_calculator(x[['MainArea[' + sample + ']'  for sample in experimental_samples_list]], \
                                       x[['MainArea[' + sample + ']'  for sample in control_samples_list]]), axis=1)
        
        # add p_value column to clean_dataframe using "p_value_calculator"
        clean_dataframe['p_val_' + self.experimental + '_' + self.control] = \
        clean_dataframe[['MainArea[' + sample + ']'  for sample in control_samples_list + experimental_samples_list]]\
        .apply(lambda x: p_value_calculator(x[['MainArea[' + sample + ']'  for sample in experimental_samples_list]], \
                                          x[['MainArea[' + sample + ']'  for sample in control_samples_list]]), axis=1)
        
        return clean_dataframe
    
    def create_volcano_plot(self, clean_dataframe):
        
        def find_fold_change_axis_limit():
            max_fold_change = np.max(clean_dataframe['fc_' + self.experimental + '_' + self.control].values)
            min_fold_change = np.min(clean_dataframe['fc_' + self.experimental + '_' + self.control].values)
            return np.ceil(max(abs(max_fold_change), abs(min_fold_change)))
        
        fold_change_axis_limit = find_fold_change_axis_limit()
        
        def find_q_value_axis_limit():
            return -np.log10(np.min(clean_dataframe['p_val_' + self.experimental + '_' + self.control].values))
            
        q_value_axis_limit = find_q_value_axis_limit()
        
        def build_volcano_df(a_lipid_class):
            # x and y values to be plotted 
            fold_change = clean_dataframe['fc_' + self.experimental + '_' + self.control][clean_dataframe['Class'] == a_lipid_class]
            pvalues = clean_dataframe['p_val_' + self.experimental + '_' + self.control][clean_dataframe['Class'] == a_lipid_class]
            
            # information to be shown when hovering over a datapoint (in addition to x and y)
            species = clean_dataframe['LipidMolec'][clean_dataframe['Class'] == a_lipid_class]
            index_list = clean_dataframe[clean_dataframe['Class'] == a_lipid_class].index
            
            return pd.DataFrame({"FC": fold_change, "qvalue": -np.log10(pvalues), "Species": species, "Index": index_list})
        
        # create a volcano plot for each lipid class 
        for lipid_class in clean_dataframe['Class'].unique():
            
            volcano_df = build_volcano_df(lipid_class)
            src = ColumnDataSource(volcano_df)
            
            # create the plot
            plot = figure(title='Volcano Plot - Class: ' + lipid_class, x_axis_label='Fold Change: ' + self.control + '/' + self.experimental,\
                          y_axis_label='q-value(-log10(p-value)')
            plot.scatter(x="FC", y="qvalue", name='volcano', size = 5, source=src)
            
            # set axis limit
            plot.x_range = Range1d(-fold_change_axis_limit, fold_change_axis_limit)
            plot.y_range = Range1d(0, q_value_axis_limit)
            
            # plot the horizontal dashed line
            plot.line(x=[i for i in range(int(-fold_change_axis_limit)-1, int(fold_change_axis_limit)+1)], y=-np.log10(0.05), line_dash = 'dashed', color='black')
            
            # plot the two vertical dashed lines
            plot.line(x=-1, y=[i for i in range(0, int(q_value_axis_limit)+1)], line_dash = 'dashed', color='black')
            plot.line(x=1, y=[i for i in range(0, int(q_value_axis_limit)+1)], line_dash = 'dashed', color='black')
                
            # add hover tool
            hover = HoverTool(tooltips = [('FC', '@FC'), ('q-value', '@qvalue'), ('Species', '@Species'), ('index', '@Index')], names=['volcano'])
            plot.add_tools(hover)
            
            # set the font size
            plot.title.text_font_size = "15pt"
            plot.xaxis.axis_label_text_font_size = "15pt"
            plot.yaxis.axis_label_text_font_size = "15pt"
            plot.xaxis.major_label_text_font_size = "15pt"
            plot.yaxis.major_label_text_font_size = "15pt"
            
            st.write('---------------------------------------------------------------------------------------------------------------------------------------------')
            st.bokeh_chart(plot)
###########################################################################################################################
st.header("LipidCruncher Demo")

st.markdown("""
        
        In the context of lipidomics, the output of mass spectrometry is the relative abundance of the lipid species that make up the sample under study. 
        This output is in the form of a spectrum in which the peaks represent an identified lipid species and the area underneath each peak reperesnts the 
        relative abundance of the corresponding lipid species. There are two pieces of software that turn this spectrum into a lipidomics dataset: 
        LipidSearch and LipidXplorer. 
        I built LipidCruncher which is a web app that allows the user to perform lipidomics analysis on the LipidSearch and LipidXplorer datasets. 
        [(Link to LipidCruncher Youtube Demo)](https://www.youtube.com/watch?v=Od9G6NYMyA0)
        
        Since the source code of LipidCruncher is few thousands lines long, in this demo, 
        I re-create a toy version of LipidCruncher with only one feature: volcano plots.  
        
        """)

st.sidebar.subheader("Upload Data")
        
lipid_search = st.sidebar.file_uploader(label='Upload your LipidSearch 4.1 dataset', type=['csv', 'txt'])
        
if lipid_search is not None:
    
    #try:
        
        df = pd.read_csv(lipid_search)
        experiment = Experiment() 
        experiment.get_input()
    
    #except:
        
        #st.sidebar.error("Pick a valid label for each condition!")
    
    #try:
        
        clean_df = CleanData(df, experiment)
        X = clean_df.clean_data()
        expand_clean_data = st.expander("View and Understand the Cleaned Dataset")
        with expand_clean_data:
            st.info("""
                    
                    Each row in the dataset represents a lipid species. There are many columns in a LipidSearch dataset, however, 
                    for the purposes of this analysis, the following columns are the ones that we need:

                    Rej: one of the built-in filtering mechanisms of LipidSearch which either takes 0 (i.e. accepted) or 1 (i.e. rejected).

                    LipidMolec: the class that the lipid species belong to and its structure (number of carbon atoms and double bonds)

                    Class: the class that the lipid species belong to

                    Calc Mass: the calculated mass of the lipid species

                    BaseRt: the retention time of the lipid species in the chromatography column

                    MainArea[s1], ..., MainArea[sN]: Area Under the Curve (AUC) representing the relative abundance of the lipid species in samples s1 to sN 
                    where N stands for the total number of the sampels
                    
                    """)
            st.write(X)
        
    #except:
        
        #st.sidebar.error("The inputs need to be modified!")
        
    #try:
    
        expand_volcano_plot = st.expander("View and Understand Volcano Plots")
        with expand_volcano_plot:
            
            st.markdown("""
                        
                        In statistics, a volcano plot is a type of scatter-plot that is used to quickly identify changes in large data sets composed of replicate data.
                        It plots significance versus fold-change on the y and x axes, respectively. 
                        A volcano plot combines a measure of statistical significance from a statistical test (e.g., a p value from a T-test) with the magnitude of the change,
                        enabling quick visual identification of those data-points that display large magnitude changes that are also statistically significant 
                        (datapoints at the top left and top right quadrant).

                        Below, q-value (i.e. -log10(p-value)) of each lipid species is plotted versus the fold change of that species. 
                        The p-value is computed from a two-sample T-test and the fold change is computed from the following formula:
                        
                        """)
                        
            latext = r'''
                    
                    $$ 
                    Fold Change = log2(\frac{Mean AUC(Condition 1)}{Mean AUC(Condition 2)})
                    $$  
            
                    '''
            
            st.write(latext)
            
            st.markdown("""
                        
                        "AUC" stands for "Area Under the Curve" which represents the relative abundance of the lipid species. 
                        "MeanAUC(condition1)" stands for the AUC averaged over all samples (replicates) belonging to condition1.
                        
                        """)
                        
            plot = VolcanoPlot(experiment)
            X = plot.add_fold_change_and_p_value_columns(X, experiment)
            plot.create_volcano_plot(X)
        
    #except:
        
        #st.error("Something went wrong!!")
            
            
