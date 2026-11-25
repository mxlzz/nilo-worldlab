# P1-D.1 — Resultado del contraste térmico acumulado
 
CI 45541091512 SUCCESS, job 138959500098, commit d4a952c23cf353566769d4d6abe3c70a500df8d0. AMD EPYC 7763, perfil heredado, script 1.864 s; no nueva mecánica ni solución FiPy, solo transporte sobre flujos archivados y control analítico. Gates vacíos, 28/28 estados, controles epoch0 B2c correctos, conservación y hashes validados. No corresponde reabrir D2 por asignación AMD: otro problema y D2 terminó inconcluso.
 
| Caso | Componente frío inicial/final (celdas, .0025 m² nominal/celda) | Déficit máximo final K | Déficit medio final K |
|---|---:|---:|---:|
| 7321_0.0 | 521/916 | 1.749970 | 0.294167 |
| 7321_0.1 | 521/921 | 3.335612 | 0.579773 |
| 8128_0.0 | 478/861 | 1.749912 | 0.279444 |
| 8128_0.1 | 478/847 | 2.759749 | 0.445516 |
 
0/28 estados con celdas significativamente más cálidas que bare (umbral -1e-7 K); residuales bare por estado <=4.657e-10 J. Frozen 7321 activo al final 850 celdas frente a 921 actual; déficit medio .527281 vs .579773 K, máximo prácticamente igual; delta campo ROI actual/frozen máximo 143.194 J por celda. Dosis total igual, pero interceptación acumulada difiere: efecto combina dosis interceptada y distribución, no aisla ventaja espacial a igual dosis sobre todo el horizonte. K0 acumula déficit por ausencia de emisión; activos mayor máximo final en parte por reducción de capacidad térmica de columna bajo erosión, no mayor refrigeración absoluta. La temperatura máxima actual está casi igual a bare (diferencia final 0.0000896/0.0003158 K activos). El área fría significa déficit relativo, no temperatura inferior a inicial. El stock granos captura energía que no llega al suelo; menor export térmico activo frente a bare (262751 vs315156 J /191454 vs227520 J), pero export material igual y no retroacción en fuerzas ni erosión. No clima calibrado, no ventaja biológica, no demostrar función frente a plantilla óptima ni novedad abierta. Archivo 14989690672, digest b580ea10b0acfdf5edb1da10118fc8faf6a17984b27155bd82af096e585416f9; gate durable pendiente. Autoverificación.
