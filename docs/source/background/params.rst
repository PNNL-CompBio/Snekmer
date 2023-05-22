.. _background-params:

Parameter Selection
===================

Snekmer's clustering mode uses various metrics, which are 
determined via user-defined parameters, in order to group 
protein sequences into clusters based on similarity. These 
clusters are thus sensitive to the parameters used to 
define similarity thresholds for sequences that are 
grouped thusly.

Why do I need to select parameters?
-----------------------------------

Different protein sequence sets will consist of proteins 
of varying lengths and similarities. Thus, an analysis
performed using one set of parameters may perform excellently
with one set of sequences, but perform poorly with another.

For instance, clustering is performed by calculating distances
or other similarity metrics between individual sequences.
Depending on the proteins involved, the magnitude of differences
between the sequences will vary, sometimes by orders of magnitude.

How should I get started with parameter selection?
--------------------------------------------------

We recommend starting with a moderate set of parameters, such
as the default parameter set provided in the template
`config.yaml <https://github.com/PNNL-CompBio/Snekmer/blob/main/resources/config.yaml>`_.
From there, depending on the clusters given by the default parameter
set, metrics such as the ``linkage`` and ``distance_threshold``, can be
further modified.