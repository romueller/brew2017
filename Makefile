##########################################################################
# This Makefile is used to control the execution of the overall analyses.
# These are started through the targets run.analysis.even respectively
# run.analysis.uneven.
#
# IMPORTANT:
# Please note that you have to provide some paths under CONFIGURATION 
# in this file as well as in cluster.sh and conf-table.sh in the scripts
# folder in order to run the analyses on your system.
# 
# The extent of the analyses (e.g. thresholds, number of replicates) can
# be changed under ANALYSIS SET-UP.
# REPLICATES is a list of replicate numbers (IDs) needing leading zeros.
# Be careful when changing the METHODS list as this might require 
# subsequent changes to cluster.sh and visualise.R in the scripts folder.
#
##########################################################################

print-%:
	@echo '$*=$($*)'

print2-%:
	@echo '$($*)'

define coloured_echo
	@tput setaf 2
	@tput bold
	@echo $1
	@tput sgr0
endef

########## CONFIGURATION ##########

# Adjust to your machine
PYTHON = python2.7 # python command / path to binary
USEARCH = ... # usearch command / path to binary


# Adjust only if you know what you are doing
FASTA_SHUFFLER = scripts/fasta_shuffler.py

###################################


######### ANALYSIS SET-UP #########

DATA_SETS = even uneven
THRESHOLDS = 1 2 3 4 5 6 7 8 9 10
REPLICATES = 0
NUM_REPLICATES = $(words $(REPLICATES))
REFS = rrna_reference.fasta

METHODS = swarm-old swarm gefast-e gefast-ef1 gefast-e2f gefast-s gefast-sf1 gefast-s2f
THREADS = 8

###################################



# ===== Data acquisition and preprocessing =====

# Format: data/<data set>.fasta
GET_DATA = $(foreach D, $(DATA_SETS), data/$D.fasta)
.SECONDEXPANSION:
$(GET_DATA): scripts/get_data.sh
	$(call coloured_echo, "Getting" $@"...")
	$(eval DATA_SET = $(subst data/,,$(basename $@)))
	bash scripts/get_data.sh $(DATA_SET)


# Format: data/<data set>.derep.fasta
DEREPLICATE_DATA = $(subst fasta,derep.fasta, $(GET_DATA))
.SECONDEXPANSION:
$(DEREPLICATE_DATA): $$(subst derep.fasta,fasta, $$@) scripts/dereplicate.sh
	$(call coloured_echo, "Dereplicating" $<"...")
	bash scripts/dereplicate.sh $< $@ 


# Format: shuffle.<data set>
SHUFFLE_DATA = $(foreach D, $(DATA_SETS), shuffle.$D)
.SECONDEXPANSION:
$(SHUFFLE_DATA): $$(subst shuffle.,, data/$$@.derep.fasta) scripts/shuffle.sh
	$(call coloured_echo, "Generating" $(NUM_REPLICATES) "shuffling(s) of" $<"...")
	bash scripts/shuffle.sh $(PYTHON) $(FASTA_SHUFFLER) $< $(NUM_REPLICATES)

# Format: data/<data set>.derep_shuffled_<replicate>.fasta
SHUFFLES = $(foreach R, $(REPLICATES), $(addsuffix _shuffled_$R.fasta, $(basename $(DEREPLICATE_DATA))))
$(SHUFFLES): $$(subst data/,shuffle., $$(basename $$(basename $$@))) ;




# ===== Taxonomic assignment =====

data/$(REFS):
	$(call coloured_echo, "Getting reference data" $<"...")
	wget -P data/ https://raw.githubusercontent.com/torognes/vsearch-eval/master/cluster/data/rrna_reference.fasta


# Format: data/<data set>.spec.ualn
ASSIGN_TAXA = $(foreach D, $(DATA_SETS), data/$D.spec.ualn)
.SECONDEXPANSION:
$(ASSIGN_TAXA): $$(subst spec.ualn,derep.fasta, $$@) data/$(REFS) scripts/taxa.sh
	$(call coloured_echo, "Assigning taxonomy for" $<"...")
	bash scripts/taxa.sh $(USEARCH) $(THREADS) $< data/$(REFS) $@




# ===== Clustering & confusion table =====

# Format: results/<data set>.metrics
METRICS_FILES = $(foreach D, $(DATA_SETS), results/$D.metrics)
.SECONDEXPANSION:
$(METRICS_FILES):
	$(call coloured_echo, "Initialising metrics file" $@"...")
	touch $@


# Requires dereplicated FASTA file (input), taxonomic assignments (input), metrics file (output)
# Format: data/<data set>.derep_shuffled_<replicate>.<threshold>.<method>.res
CLUSTER_CONF = $(foreach T, $(THRESHOLDS), $(foreach M, $(METHODS), $(addsuffix .$T.$M.res, $(basename $(SHUFFLES)))))
.SECONDEXPANSION:
$(CLUSTER_CONF): $$(addsuffix .fasta, $$(basename $$(basename $$(basename $$@)))) \
		$$(addsuffix .spec.ualn, $$(basename $$(basename $$(basename $$(basename $$@))))) \
		$$(addsuffix .metrics, $$(subst data/,results/, $$(basename $$(basename $$(basename $$(basename $$@)))))) \
		scripts/cluster.sh scripts/conf-table.sh
	$(eval FASTA = $<)
	$(eval TAXA = $(addsuffix .spec.ualn, $(basename $(basename $(basename $(basename $@))))))
	$(eval THRESHOLD = $(subst .,, $(suffix $(basename $(basename $@)))))
	$(eval METHOD = $(subst .,, $(suffix $(basename $@))))
	$(eval SHUFFLING = $(subst .derep_shuffled_,, $(suffix $(basename $(basename $(basename $@))))))
	$(eval CONF = $(addsuffix .conf, $(basename $@)))
	$(eval METRICS = $(addsuffix .metrics, $(subst data/,results/, $(basename $(basename $(FASTA))))))
	$(call coloured_echo, "Clustering" $(FASTA) "with" $(METHOD) "(replicate" $(SHUFFLING)", threshold" $(THRESHOLD)")...")
	bash scripts/cluster.sh $(FASTA) $@ $(THRESHOLD) $(TAXA) $(METHOD) $(THREADS)
	bash scripts/conf-table.sh $@ $(SHUFFLING) $(METHOD) $(THRESHOLD) $(TAXA) $(CONF) $(METRICS)
	rm $@ $(CONF)




# ===== Plotting =====

# Format: results/<data set>.metrics.long
PRODUCE_LONG_TABLE = $(addsuffix .long, $(METRICS_FILES))
.SECONDEXPANSION:
$(PRODUCE_LONG_TABLE): $$(basename $$@) scripts/produce_long-form.sh
	$(call coloured_echo, "Producing long table from" $<"...")
	bash scripts/produce_long-form.sh $< $@

COMMA :=,
SPACE :=
SPACE +=

# Note: All plot functions use the long form of the metrics table.

# Format (of all plot targets): results/<data set>.<plot type>.pdf
PLOT_METRICS = $(addsuffix .pdf, $(METRICS_FILES))
.SECONDEXPANSION:
$(PLOT_METRICS): $$(subst pdf,long, $$@) scripts/visualise.R
	$(call coloured_echo, "Plotting " $@"...")
	$(eval METRICS_FILE = $<)
	$(eval THRESHOLD_LIST = $(subst $(SPACE),$(COMMA),$(strip $(THRESHOLDS))))
	R -e "source('scripts/visualise.R'); plot_metrics('$(METRICS_FILE)', '$@', c($(THRESHOLD_LIST)))"

PLOT_MEDIANS = $(subst .metrics,.medians.pdf, $(METRICS_FILES))
.SECONDEXPANSION:
$(PLOT_MEDIANS): $$(subst medians.pdf,metrics.long, $$@) scripts/visualise.R
	$(call coloured_echo, "Plotting " $@"...")
	$(eval METRICS_FILE = $<)
	$(eval THRESHOLD_LIST = $(subst $(SPACE),$(COMMA),$(strip $(THRESHOLDS))))
	R -e "source('scripts/visualise.R'); plot_medians('$(METRICS_FILE)', '$@', c($(THRESHOLD_LIST)))"

PLOT_MEDIANS_JITTER = $(subst .metrics,.medians-jitter.pdf, $(METRICS_FILES))
.SECONDEXPANSION:
$(PLOT_MEDIANS_JITTER): $$(subst medians-jitter.pdf,metrics.long, $$@) scripts/visualise.R
	$(call coloured_echo, "Plotting " $@"...")
	$(eval METRICS_FILE = $<)
	$(eval THRESHOLD_LIST = $(subst $(SPACE),$(COMMA),$(strip $(THRESHOLDS))))
	R -e "source('scripts/visualise.R'); plot_medians_jitter('$(METRICS_FILE)', '$@', c($(THRESHOLD_LIST)))"

PLOT_MEDIANS_LINE = $(subst .metrics,.medians-line.pdf, $(METRICS_FILES))
.SECONDEXPANSION:
$(PLOT_MEDIANS_LINE): $$(subst medians-line.pdf,metrics.long, $$@) scripts/visualise.R
	$(call coloured_echo, "Plotting " $@"...")
	$(eval METRICS_FILE = $<)
	$(eval THRESHOLD_LIST = $(subst $(SPACE),$(COMMA),$(strip $(THRESHOLDS))))
	R -e "source('scripts/visualise.R'); plot_medians_line('$(METRICS_FILE)', '$@', c($(THRESHOLD_LIST)))"

PLOT_SINGLE_METRIC_JITTER = $(subst .metrics,.single-metric-jitter.pdf, $(METRICS_FILES))
.SECONDEXPANSION:
$(PLOT_SINGLE_METRIC_JITTER): $$(subst single-metric-jitter.pdf,metrics.long, $$@) scripts/visualise.R
	$(call coloured_echo, "Plotting " $@"...")
	$(eval METRICS_FILE = $<)
	$(eval THRESHOLD_LIST = $(subst $(SPACE),$(COMMA),$(strip $(THRESHOLDS))))
	R -e "source('scripts/visualise.R'); plot_single_metric_jitter('$(METRICS_FILE)', '$@', c($(THRESHOLD_LIST)))"

PLOT_SINGLE_METRIC_LINE = $(subst .metrics,.single-metric-line.pdf, $(METRICS_FILES))
.SECONDEXPANSION:
$(PLOT_SINGLE_METRIC_LINE): $$(subst single-metric-line.pdf,metrics.long, $$@) scripts/visualise.R
	$(call coloured_echo, "Plotting " $@"...")
	$(eval METRICS_FILE = $<)
	$(eval THRESHOLD_LIST = $(subst $(SPACE),$(COMMA),$(strip $(THRESHOLDS))))
	R -e "source('scripts/visualise.R'); plot_single_metric_line('$(METRICS_FILE)', '$@', c($(THRESHOLD_LIST)))"




# ===== Overall analyses =====

run.analysis.even: $(filter-out data/uneven%, $(CLUSTER_CONF)) \
			$(filter-out results/uneven%, $(PLOT_METRICS)) \
			$(filter-out results/uneven%, $(PLOT_MEDIANS)) \
			$(filter-out results/uneven%, $(PLOT_MEDIANS_JITTER)) \
			$(filter-out results/uneven%, $(PLOT_MEDIANS_LINE)) \
			$(filter-out results/uneven%, $(PLOT_SINGLE_METRIC_JITTER)) \
			$(filter-out results/uneven%, $(PLOT_SINGLE_METRIC_LINE))


run.analysis.uneven: $(filter data/uneven%, $(CLUSTER_CONF)) \
			$(filter results/uneven%, $(PLOT_METRICS)) \
			$(filter results/uneven%, $(PLOT_MEDIANS)) \
			$(filter results/uneven%, $(PLOT_MEDIANS_JITTER)) \
			$(filter results/uneven%, $(PLOT_MEDIANS_LINE)) \
			$(filter results/uneven%, $(PLOT_SINGLE_METRIC_JITTER)) \
			$(filter results/uneven%, $(PLOT_SINGLE_METRIC_LINE))

