The data includes PTID of 635 subjects, thickness measured at baseline, 5 covariates and 8 SNPs, and disease status (AD/MCI/CN).

1) column1 : subject ID

2) column2 - column69 : 68 ROIs (Desikan-Killany atlas (Desikan et al.2006))
   column2 - column13 : ROIs for default mode network
		"ST31TA", "ST90TA",  # Inferior parietal
		"ST32TA", "ST91TA",  # Inferior temporal
		"ST39TA", "ST98TA",  # Medial orbitofrontal (medial prefrontal cortex)
		"ST44TA", "ST103TA", # Parahippocampal
		"ST52TA", "ST111TA", # Precuneus
		"ST50TA", "ST109TA", # Posterior cingulate

   column13 - column69 : rest of ROIs
		"ST13TA", "ST72TA",  # Banks of the superior temporal sulcus
		"ST14TA", "ST73TA",  # Caudal anterior cingulate
		"ST15TA", "ST74TA",  # Caudal middle frontal
		"ST23TA", "ST82TA",  # Cuneus
		"ST24TA", "ST83TA",  # Entorhinal
		"ST25TA", "ST84TA",  # Frontal pole
		"ST26TA", "ST85TA",  # Fusiform
		"ST34TA", "ST93TA",  # Isthmus of the cingulate
		"ST35TA", "ST94TA",  # Lateral occipital
		"ST36TA", "ST95TA",  # Lateral orbitofrontal
		"ST38TA", "ST97TA",  # Lingual
		"ST40TA", "ST99TA",  # Middle temporal
		"ST43TA", "ST102TA", # Paracentral
		"ST45TA", "ST104TA", # Pars opercularis
		"ST46TA", "ST105TA", # Pars orbitalis
		"ST47TA", "ST106TA", # Pars triangularis
		"ST48TA", "ST107TA", # Pericalcarine
		"ST49TA", "ST108TA", # Postcentral
		"ST51TA", "ST110TA", # Precentral
		"ST54TA", "ST113TA", # Rostral anterior cingulate
		"ST55TA", "ST114TA", # Rostral middle frontal
		"ST56TA", "ST115TA", # Superior frontal
		"ST57TA", "ST116TA", # Superior parietal
		"ST58TA", "ST117TA", # Superior temporal
		"ST59TA", "ST118TA", # Supramarginal
		"ST60TA", "ST119TA", # Temporal pole
		"ST62TA", "ST121TA", # Transverse temporal
		"ST129TA","ST130TA"  # Insula
			

3) column70 - column75 : 5 covariates (PTGENDER, AGE, PTEDUCAT, PTHAND, ICV_bl)

4) column76 - column82 : genotype for 8 SNPs (used in Shen et al.2010 Table4)
   APOE4, rs2075650, rs7526034, rs10932886, rs7647307, rs7610017, rs4692256, rs6463843

5) column83: disease status: AD, MCI, or CN (normal controls)

	

