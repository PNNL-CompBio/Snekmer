Introduction to Machine Learning
================================

Snekmer uses the AAR-kmer approach to develop machine learning 
models for the characterization of protein sequences. To
familiarize users with the terminology and ideas underlying
Snekmer, we will include a brief introduction to machine learning.

What is Machine Learning?
-------------------------

Machine learning (ML) encompasses many methods by which computers,
or "machines," are taught to recognize, or "learn," underlying
patterns in data. Essentially, a machine learning algorithm is
trained on a given dataset for a particular task. Machine learning
has been employed to solve a diverse array of problems, including
detecting spam messages,\ :footcite:p:`GmailSpam`
translating text,\ :footcite:p:`NLLBTeam2022` creating
images,\ :footcite:p:`Ramesh2022`
and more.

Depending on the type of problem, machine learning can take a few
different forms. The primary categories of machine learning are
*supervised learning*, *semi-supervised learning*, and 
*unsupervised learning*.

Supervised Learning
:::::::::::::::::::

In **supervised learning**, the user supplies a dataset with known
labels, and then trains an algorithm to accurately predict the
labels. The resulting model is deemed *supervised*, since the model
learns relationships that have already been determined by humans,
as though humans are "supervising" the learning task.

For instance, Snekmer's ``snekmer model`` mode uses supervised
machine learning. The models produced by Snekmer are created by
first quantifying how AAR-kmers are distributed across pre-annotated
protein groupings and subsequently training a machine learning
model to learn the relationship between AAR-kmer distributions
and a given annotation.

Unsupervised Learning
::::::::::::::::::::::

In **unsupervised learning**, the user supplies a dataset where the
underlying labels or relationships are not known. Thus, the
resulting model is trained to identify relationships or patterns in
the data independent of known labels. These models are called
*unsupervised* in contrast with supervised machine learning models;
the algorithm must determine its own assessments of relationships
in the underlying data without labels that have been pre-drawn
by humans.

Snekmer's ``snekmer cluster`` mode is an example of an unsupervised
machine learning. In clustering mode, Snekmer trains clustering
models to group protein sequences based on their similarity to other
sequences in the dataset in accordance with AAR-kmer profiles.
These groupings can be based on whichever distance or similarity
metrics are specified by the user. Thus, the parameters used to
examine the protein sequences are critical in calculating clusters.
Changing the AAR alphabet, kmer length, or clustering method can
have a large effect on how the algorithm assesses sequences.

Semi-Supervised Learning
::::::::::::::::::::::::

As its name may suggest, **semi-supervised learning** lies in
between supervised and unsupervised machine learning. In semi-
supervised learning, the user supplies a small corpus of
labeled data alongside a larger set of unlabeled data. The machine
learning model is trained on both labeled and unlabeled data points.
Semi-supervised learning is useful in cases where producing a large
dataset for supervised machine learning is not realistic.

Citations
:::::::::

.. footbibliography::