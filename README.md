Implementing realistic grazing in models of catastrophic shift (CAS01)
======================================================================

## Objective

This project aims to provide more realistic models for the resilience of dryland systems. 
In comparison to the parent facilitation model used by Kéfi et al. 2007, two additive assumptions were explored:  

 - 	**associative protection** assumes that grazers affect isolated plants more than plants growing associatively, in patchy structures. This adds a second facilitative feed-back to the parent model of Kéfi et al. on the mortality of the plant individual. At high vegetation cover, less plants are vulnerable to grazing. As a consequence, the vegetation structure is altered. The unit of a patch becomes consolidated. Hypothesis: the patchiness of the landscape increases. the vegetated state becomes more consolidated.  
 - 	**livestock model** assumes that the feeding on vegetation depends on the vegetation cover (functional response) and thus the number of plants dying due to grazing increases with vegetation cover. While the parent model assumes a constant risk for each individual plant, this model renders the individual risk dynamically: It is highest for very low cover and approaches zero when cover is high.  

One additional dimension is added in Alain's project:  

 - 	**palatability model** assumes that plant species differ in vulnerability to grazing which would cause very different stable states.  

### Pressure
The two models mimic well founded mechanisms in pastural landscapes and can be applied in cellular automata models. The model simulations give a clue about how these mechanisms change the dependency of the landscapes resilience on two types of pressure: environmental quality and grazing intensity. 

- **environmental quality** boils down to the most relevant limiting factor in drylands: Water availability. If the landscape suffers reduced precipitation or higher temperatures, the vegetation responds nonlinearly. Facilitation improves environmental quality locally, and thus provides a characteristic spatial vegetation structure. This structure collapses as the environment becomes more arid. This is the major driver behind the bistability of the dryland model. In the model, local facilitation improves the regeneration of degraded cells and the recolonisation of empty cells. 

-	**grazing intensity** is increasing plant's mortality. Different definitions of grazing intensity are imaginable: It can be a constant effect on the individual plant's mortality (parent model) or can be a pressure applied to the entire landscape (livestock) that is distributed among all plant individuals. It can be formulated to affect all plants equally, or to affect plants differently depending on the local neighborhood (associative protection). While the model so far assumes all plants equal, grazers would prefer some species over others, which produces very different spatial vegetation structure. Finally, on the patch level, grazing can be formulated to affect all patches equally or to include a preference of grazers for large patches to reduce search times.  
	In reality, a number of grazer individuals move through the landscape and make choices on the patch scale as well as on the individual plant scale (Adler 20xx). It would graze a certain amount of vegetation each day, while preferring easily accessible plants. Thus, an individual grazer model would include all of the previous assumptions.
 
Both types of pressure affect different life history stages of the plant individual and additively determine the formation of patches. 

### Landscape response
The numerical simulation of the different models, over gradients of the two types of pressure allow investigations of the system resilience and of indicators of catastrophic shifts.  

- **system resilience** can be looked at in two ways: The overlap of the alternative stable states (vegetated vs. desert) and the robustness of the stable states to perturbations (relation between stable and unstable equilibria).

- **indicators** can be distinguished in spatial indicators, looking at patch patterns on the landscape scale, and temporal indicators, looking at the change and variability of metrics over time. The spatial indicators include patchiness, size of largest patch, as well as descriptions of cumulative patch size distributions. 

### Questions 

**How does increasing grazing 'realism' affect ecosystem resilience?**


**How do different types of local facilitation affect ecosystem resilience and spatial patterns?**

Biotic vs. abiotic facilitation mechanisms. Both affect the emergent patch structure, through different parts of plant life history. The **associative protection model** seems suitable to investigate how the system responds to stress in terms of resilience and spatial structure. 

**How do spatial explicit pressure affect the applicability of early warning signs?**

Ecosystems under stress respond more or less quickly with changes in spatial structure. 


### Deliveries

We think about splitting the project in the following deliveries.
 
#### Manuscript: "Grazing raises risk of catastrophic shifts"

#### Manuscript: "Biotic vs abiotic facilitation as drivers of ecosystem resilience and spatial structure"



## License

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">Towards a spatial resilience framework for arid and semi-arid rangelands</span> by <a xmlns:cc="http://creativecommons.org/ns#" href="https://fdschneider.github.io" property="cc:attributionName" rel="cc:attributionURL">Florian D. Schneider</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/cascade-wp6/resilience_rangelands" rel="dct:source">https://github.com/cascade-wp6/resilience_rangelands</a>.
