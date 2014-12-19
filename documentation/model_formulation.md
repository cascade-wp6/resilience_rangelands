# Basic Model formulation
Florian D. Schneider  
Tuesday, November 18, 2014  

This document is supposed to become an interactive model representation of the pair-approximation model. 

The default form of the differential equation describing vegetation cover is 

$$ \frac{dV}{dt} = G(V) - C(V)H $$

with growth function $G$ and a function $C$ which represents consumption by a herbivore multiplied by the density of herbivores, $H$. Formally, intrinsic mortality could be included in both terms. In the formulations of @noy-meir75 and @vandekoppel97 it is assumed to be covered by the logistic growth function.

It might become important in the pair-approximation model, where intrinsic plant mortality is included as a proportional loss of mature plants.



## non-spatial models

### logistic growth

If assuming phenomenological logistic growth, the gain in vegetation as a function of present vegetation cover is 

$$ G(V) = rV(1-V/K) \,,$$

with the growth rate, $r$, and the carrying capacity, $K$, decelerating the exponential growth as vegetation cover increases and reducing it to zero if $V \geq K$. Thus, it is a decelerating effect (a negative feed-back) at high vegetation cover. 

![plot of chunk unnamed-chunk-2](model_formulation_files/figure-html/unnamed-chunk-2.png) 

#### mechanistic growth assumptions

The logistic growth model is phenomenological and does not provide a particular ecological explanation, why growth is decelerating to zero with increasing cover. It assumes just some form of global competition leading to a balance of gains and losses at the carrying capacity of the ecosystem. 

Mechanistic explanations can be the limitations by a globally diffuse resource, like water, or the occupation of space and the subsequent limitation of photosynthetic active surface. 

A multitude of such global feedback processes can define the precise shape of the growth function leading to asymetrical or truncated shapes, which could importantly affect the bistability properties of the system [@noy-meir75]. He defines

> The general shape of the growth function relating net growth to plant biomass, $G(V)$, is convex with a single maximum (Fig. 1). The increase in the low-biomass range expresses the increase in photosynthetic capacity with increasing leaf area. The levelling
off to a plateau and the decrease in the higher range expresses self-interference effects (shading, competition); at the point $V = V_\mathrm{max}$ 'maintenance' losses equal photosynthesis and $G = 0$. Vmax is the maximum (stable) biomass in ungrazed vegetation. [...] Of the explicit forms the best known is the logistic [...] but other functions are possible.


![plot of chunk unnamed-chunk-3](model_formulation_files/figure-html/unnamed-chunk-3.png) 

These are global feed-backs. But local feedbacks also could lead to this deceleration:  
For instance, if plant cover would have a dominantly negative effect on it's surroundings, by depleting nutrients or light, the local neighborhood of a plant cannot be occupied as easy as more remote places. As cover increases and the interspace areas are closing in, the space available for rejuvenation approaches zero. 

#### decelerating growth at low cover  

In arid ecosystems, growth would be decelerated at very low cover, because retention of water and organic matter in the site is fundamentally reduced. This postive feed-back (low cover reduces growth, leading to even lower cover) was discussed by @vandekoppel97 and @rietkerk97 to explain bistability in arid ecosystems.

$$ G(V) = rV(1-V/K) \,,$$


### functional response

#### linear plant mortality

In the cellular automata model and mean-field approximations of @kefi07, individual plant mortality is a constant rate, which means it is increasing linearly with vegetation cover. 

Plant mortality due to grazing therefore would be described as

$$ C(V) = mV \,, $$

with the grazing mortality, $m$.

![plot of chunk unnamed-chunk-4](model_formulation_files/figure-html/unnamed-chunk-4.png) 

#### constant plant mortality

In fact, when assuming a constant density of herbivores, grazing mortality will certainly saturate at some point, due to limitations of the herbivores handling times. A very parsimonious model would assume that grazing always is defined by this limitation of handling times, even at very low densities. 

Mortality due to grazing,

$$ C(V) = m \,, $$

thus would be defined as constant, $m$, regardless of the vegetation cover.


![plot of chunk unnamed-chunk-5](model_formulation_files/figure-html/unnamed-chunk-5.png) 


#### type 2 functional response

In natural systems, however, consumption is defined by saturating functional responses, i.e. at low cover, density dependent mortality defined by search efficiency dominates, but at high cover handling times of the consumer will limit its intake rates. 

Therefore, plant mortality due to grazing should be 

$$ C(V) = \frac{aV}{1+ahV} \,,$$

with the handling time, $h$ [time per resource unit], and the search efficiency, $a$ [area foraged in per time].  

![plot of chunk unnamed-chunk-6](model_formulation_files/figure-html/unnamed-chunk-6.png) 


#### type 3 functional response

$$ C(V) = \frac{bV^{1+q}}{1+bhV^{1+q}} $$


![plot of chunk unnamed-chunk-7](model_formulation_files/figure-html/unnamed-chunk-7.png) 



### bistability

Plotting the growth rate against the different mortality rates shows how the definition of grazing mortality alters our prediction of the steady states of the system. 

Compared to the more realistic non-linear functions, the constant grazing mortality assumption seems to predict both the attractor of the vegetated state and the tipping point much better than the density dependent mortality used in Kéfi et al 2007 TPB. In fact, the linear density dependent mortality does not predict bistability at all, given a simple logistic growht function. 

The critical threshold that would lead to degradation is overestimated by the constant model. 

The type 3 functional response finally predicts the existance of a stable degraded state that is not a desert, i.e. not entirely unvegetated. 



![plot of chunk unnamed-chunk-8](model_formulation_files/figure-html/unnamed-chunk-8.png) 


### open questions

#### Intrinsic mortality

The mortality of plants regardless of grazing cover must be considered to be already internalised into the logistic growth definition. 

Alternatively, it can be considered as an independent constant rate, $M(V) = mV$, 

$$ \frac{dV}{dt} = G(V) - M(V) - C(V) $$

and thus would affect the alternative stable states of the ecosystem.


![plot of chunk unnamed-chunk-9](model_formulation_files/figure-html/unnamed-chunk-9.png) 



## spatially-explicit models

The models described above were discussed by @noy-meir75 for grazing systems where vegetation is assumed to be of homogeneous distribution, and no spatially-explicit effects are of relevance. 
Refugia or fencing of areas, however, could be considered as spatial effects. 

In arid and semi-arid rangelands, spatial vegetation structure was considered to be of major importance for plant growth and mortality [@rietkerk97; @kefi07]. It generates patterns that act as positive feed-backs on the local scale by providing refugia to grazing and overcome growth limitations at low vegetation cover.

These mechanisms can be implemented in the model by using a pair-approximation approach. 

### growth

In the formalism of the pair-approximation model of @kefi07, growth, *i.e.* the total amount of empty cells transitioning into vegetated cells, becomes a complex function of the relative cover of the multiple pairs including a $+$.

growth part of equation (28) from @kefi07:

$$  \rho_{0}w_{0,+} = \left( \delta \rho_+ + (1-\delta) \frac{\rho_+ - \rho_{+-} - \rho_{++} }{1 - \rho_+ - \rho_-} \right) (b - c \rho_+)(1- \rho_+ - \rho_-)  $$

The goal of this exercise would be to transform this equation into a function of $\rho_+$, assuming that the landscape is at a vegetated equilibrium. 

Therefore, we are looking for the general case of $G(\rho_+)$ where

$$ \frac{d\rho_+}{dt} = \frac{d\rho_-}{dt} = \frac{d\rho_{+,-}}{dt}  = \frac{d\rho_{+,+}}{dt} = \frac{d\rho_{-,-}}{dt}   = 0$$

Equilibrium has to fulfil

$$ \rho_0w_{0,+} = - m \rho_+ \,.$$


thus

$$  - m \rho_+ = \left( \delta \rho_+ + (1-\delta) \frac{\rho_+ - \rho_{+-} - \rho_{++} }{1 - \rho_+ - \rho_-} \right) (b - c \rho_+)(1- \rho_+ - \rho_-)  $$



### spatially-explicit grazing mortality




# References
