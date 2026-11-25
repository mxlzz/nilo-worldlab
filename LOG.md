# Bitácora

- 2026-11-24: usuario aprobó diseño y migración. Origen inspeccionado, destino inicialmente vacío; corte 0448ea6e5e02d77142cec270ea99ac0773493fc1.
- Se recuperaron 12 commits desde dc1fdff. P0-D/P0-E completos, P1 aprobada, P1-A original negativa instrumental; P1-A-R CI 45488606315 y log de 12 casos con gates satisfechos. No se ejecutó física nueva. P1-B es siguiente fase científica.
- Se inició traspaso con workflow exclusivo de copia, comprobación de blobs y archivo de evidencia. Estado requiere lectura verificada y cobertura antes de avanzar.
- WL-M2: transferencia 45491750969 SUCCESS; 314/314 blobs exactos, 58/58 artefactos digest verificado. Origen permanece en 0448ea6e5e02d77142cec270ea99ac0773493fc1.
- WL-M3: archivo de logs 45492349703 SUCCESS, 481/487 ZIP, seis 404 sin jobs. Manifiestos y archivos recuperados en evidence/. Ningún artefacto inventariado faltante.
- WL-M4: verificación de árboles, evidencia, documentos, enlaces y simulación de continuidad (evidence/MIGRATION_VERIFICATION.md). Se archivan workflows de traspaso. P1-B.1 es siguiente unidad científica; STATE se actualizará al final.
- P1-B.1: inspección de contratos y fuentes Landlab v2.11.0; identificado déficit de flujos brutos en H33. Prerregistro de receptor por partición del suelo existente, entalpía incremental y transporte sedimentario. GO documental condicionado a B.2a; sin ejecución física nueva. Ver research/P1-B/P1-B.1_CONTRATO_Y_PRERREGISTRO.md.

- P1-B.2a: implementado adaptador de solo lectura de campos nativos SPACE (tasas × dt × área/densidad según unidades), comparador de hashes P0-D, controles y CI de 15 min. Referencias archivadas con digest verificado al ejecutar. Lanzamiento por push; resultado pendiente, no habilitar B.2b.

- P1-B.2a CI 45502145469 SUCCESS, gates vacíos, 28/28 hashes exactos y 72 pasos, autoverificación. Se inicia archivo durable del artefacto y log antes de B.2b.

- Archivo B.2a run 45503207288: ZIP pasó digest; gh CLI rechazó secuencias ANSI de log (no se escribió evidencia binaria). Recuperación limitada: habilitar salida cruda exclusivamente al archivo de log, sin imprimirlo; checkout con historial para rebase. No repetir física.

- Archivo recuperado CI 45503283036 SUCCESS: ZIP, log, JSON y manifest en evidence/P1-B/B2a, digest y NPZ verificados en CI, tamaños y blobs releídos en HEAD 717caa02. B.2a cerrada, workflow temporal de archivo inactivo. GO B.2b.

- P1-B.2b: operador de energía de partición disjunta, casos analíticos, rechazo corrupción, remapeo, contraste por nodo con ConcentrationTrackerForSpace. Owner reducido de laboratorio, no integración P0-D ni dinámica térmica. Se usa entalpía específica inicial uniforme 1 J/kg y gradiente determinista con ROI heterogénea (solo si pasa uniforme); 18 pasos en cuatro casos archivados. CI limitado a 15 min.

- P1-B.2b CI 45503529832 SUCCESS, gates vacíos, ocho casos 144 pasos, tracker max 4.191e-9 J/nodo, residual 2.3842e-7 J, propietario reducido. Archivo de evidencia pendiente antes de B.2c.

- Archivo B.2b CI 45503722111 SUCCESS, ZIP/log/JSON/manifest verificados, blobs y tamaños releídos en b87b19cf. Workflow de archivo desactivado. B.2b cerrada en esquema reducido, GO B.2c.

- P1-B.2c implementada según gates prerregistrados: kernel B.2b copiado literalmente y verificado al iniciar, cuatro propietarios P0-D, rayos exactos/FiPy, control/nulos, transporte de tres pasos archivados, remapeo y resolución, owner completo de snapshot. Sin repetir mecánica ni Landlab. CI pendiente; no declarar resultado hasta lectura.

- P1-B.2c CI 45504052684 SUCCESS, 4 casos, gates vacíos, 521/478 celdas vs P95 9/12; energía exportada y retenida contabilizada. Cierre P1-B condicionado al archivo durable. No afirmar dinámica viva ni E7.

- Archivo B.2c CI 45516061554 SUCCESS, ZIP/log/JSON/manifest verificados, los 44 outputs cotejados por digest y blobs/tamaños releídos en 62ca60d5. P1-B cerrada limitada. P1-C.1 prerregistrada; C.2a siguiente unidad.

- C.2a implementado productor continuo: worker por demanda con gates B.2a, owner P1C-v1, irradiación acumulada, controles observador/flujo/pose y checkpoint epoch3. Solo referencia posterior a evolución, no modificar inherited. Lanzamiento CI 15 min; aún sin resultado.

- C.2a intento 45516385783 FAILURE: conteo incorrecto 12 vs 11 campos, confirmado en 4 casos. Archivar antes de único reintento. No gates físicos modificados.

- C.2a fallo archivado CI 45516870053 SUCCESS, ZIP/log/JSON y manifest verificados; corregido trigger de archivo sin repetir física. Revisión de comparador: reemplazo de conteo erróneo por conjunto exacto de 33 nombres para tres pasos, mantiene comparación numérica y topológica. Único reintento instrumental C.2a, sin cambiar parámetros ni criterios físicos.

- C.2a segundo intento CI 45516933103 SUCCESS; cuatro casos y tres controles, 28 proyecciones, 72 pasos, exactitud y balances; se inicia archivo separado producer/checkpoint y log antes de C.2b.

- Archivo C.2a CI 45518056825 SUCCESS, ambos ZIP, 147 outputs, 12 archivos de checkpoint, log/JSON idénticos verificados; blobs releídos en 3918169f. C.2a aprobada, GO C.2b sin entregar futuro al consumidor. Workflow temporal inactivo.

- C.2b implementada: consumer sparse checkout código/raw/checkpoint únicamente, .git removido y lista de evidencia exacta, sin datos futuros del productor. Reconstrucción nativa P0-D sin pasos previos, siete controles de corrupción, 9 soil/2160 Bullet/18 térmicos nuevos por caso. Comparador separado posterior verifica hashes lógicos, arrays bit a bit y reporta igualdad de bytes. CI 15+5 min, pendiente.

- C.2b CI 45518264979: consumidor pasa 4 casos/28 corrupciones; comparador falla exactitud owner 12/12 pese a trayectorias/flujos/rayos exactos. Archivo de ambos ZIP/logs y diagnóstico no mutante, sin nueva dinámica. NO-GO P1-D.

- Archivo C.2b intento 45519115487 sin jobs, YAML con :all: en escalar plano; corregido a bloque shell. Solo archivo/diagnóstico, no repetición científica.

- Archivo C.2b CI 45519166290 SUCCESS, ambos ZIP/log/JSON, diagnóstico 12 estados confirmado; diferencias reales E_roi hasta 9.095e-13 J, E_rest hasta 7.276e-12 J. NO-GO P1-D, prerregistro D2 mínimo sin dinámica. Archivo temporal desactivado.

- D2 implementado y lanzado: 2 jobs x 2 perfiles x 2 repeticiones en dos casos retenidos; un solo bloque térmico, sin dinámica nueva, perfil candidato Nehalem predefinido. Gate C.2b intacto.

- D2 CI 45519386325 SUCCESS técnico pero inconcluso científico: dos Intel, no AMD; native reproduce referencia, Nehalem cambia hashes reproduciblemente sin cambiar transporte. Agotado presupuesto D2; no habilitar C.2b-R. Archivo pendiente; elevar decisión de alcance.

- D2 archivo CI 45520590725 SUCCESS, ZIP/log/JSON y 16 outputs NC cotejados, blobs releídos en fff200ab. Dos runners Intel, sin AMD; API de runners no accesible desde conector, disponibilidad dirigida no confirmada. Detener en decisión D2: preservar exactitud y autorizar verificación AMD dirigida/acotada, o redefinir alcance explícitamente, o cerrar P1 parcial. No nuevo lanzamiento aleatorio ni P1-D.

- Usuario autorizó extensión única D2-AMD sin gasto. Runner público estándar sin selección de fabricante documentada; una asignación con guard AMD antes de instalar, no reintentos si Intel. Protocolo D2 intacto.

- D2-AMD CI 45521624938: Intel, filtro falló antes de ejecución científica, no AMD. Única asignación consumida, no relanzar. Archivar CPU/log y elevar disponibilidad/alcance.

- D2-AMD archivo CI 45521720220 SUCCESS, CPU/log/ZIP verificados, blobs releídos en 7a8113a9. Desactivados workflows de asignación única y archivo. Bloqueo disponibilidad AMD, no fallo científico nuevo. Requiere acceso a AMD dirigido o decisión explícita de alcance; no reintentar asignaciones aleatorias.

- Usuario aprueba separar persistencia exacta en entorno especificado y portabilidad numérica. Prerregistro S1/S2: misma referencia, mismo consumidor, CPU Intel equivalente; S2 evaluación retrospectiva explícita de AMD, límites físicos predefinidos, sin generalización. No aprobar C.2b original retrospectivamente.

- S1/S2 lanzados: consumidor y comparador originales intactos, guard CPU Intel, portabilidad retrospectiva solo después del comparador exacto, controles 1J, L1 sin cancelación. Sin repetir productor ni modificar leyes.

- S1 CI 45531989960 inconclusa: Intel 8573C distinto, guard sin física; S2 omitida. Log y jobs archivados. Prerregistro S1-H: referencia nueva y procesos separados mismo host, sin selección CPU; gates intactos, comparación posterior, S2 histórica. No repetir asignaciones aleatorias.

- S1-H lanzada: original productor/consumer/compare sin modificar; nueva referencia y checkpoint empaquetados y publicados, futuro local eliminado, consumidor separado, comparación en job posterior y S2 histórica. Sin selección CPU ni ajustes de gates.

- S1-H/S2 CI 45539671147 SUCCESS, 12/12 sufijos exactos, 36/36 archivos byte-exactos, 28 corrupciones, S2 12/12 y 24 controles; máximo L1 1.000102e-10 J. Archivo durable pendiente antes de P1-D.

- S1-H archivo CI 45540857670 SUCCESS, 5 ZIP y 3 logs, 147 producer/12 checkpoint/36 suffix hashes cotejados; blobs/tamaños releídos en 18b4d446. Nueva referencia coincide además 147/147 archivos con productor original. P1-C cerrada limitada según alcance aprobado, C2b original negativa, S1 original inconclusa, S2 retrospectiva; GO P1-D.

- P1-D.1 prerregistro de auditoría: contraste térmico de 28 snapshots vivos frente a bare analítico uniforme, mismos flujos y stock, controles B2c y frozen retenido, sin repetir mecánica ni solver térmico. Evaluación descriptiva y gate de instrumento.

- P1-D.1 lanzado: lectura de owners y flujos por digest, bare analítico uniforme con operador conservativo heredado, contraste por época y frozen. Export material cotejado con tolerancia heredada .09216 kg de C.1, sin cambios. CI 5 min, no nuevos solvers ni mecánica.

- P1-D.1 CI 45541091512 SUCCESS, 28 contrastes, 4 casos, 7 frozen. Déficit final activo 3.3356/2.7597 K máximo, 921/847 celdas conectadas, controles B2c y balance superados. Script 1.864s, archivo pendiente. Sin inferencia E7 ni retroacción bidireccional.

- P1-D.1 archivo CI 45541914771 SUCCESS, ZIP/log/JSON digest y blobs verificados. P1-D.2 dictamen: P1 cerrada como consecuencia ambiental causal finita; no retroacción bidireccional ni E7 general. GO diseño E7.1 prospectivo discriminante, no repetición de montículos por inercia, no E8. Encargo y límites registrados.

- E7.1 diseño: NO-GO nuevas sombras; GO historia de suministro con idénticos puntos y distribución de reacción basal como función distinta. Antecedentes Vanel/Atman; reutilizar Bullet/H35 y terreno retenido. E7.1a calibrar fuerza completa y CDF antes de formación prospectiva; E7.1b aún no prerregistrada ni ejecutada.

- E7.1a instrumento implementado y lanzado: cuatro fixtures de reacción basal, sumatoria normal+tangencial, CDF sintético y poses no mutadas, 240 muestras por fixture; sin escenas prospectivas. CI hasta5min gratuito.

- E7.1a CI 45595575928 SUCCESS: 4/4 fixtures, CDF sintético, reacción y fricción inclinada, residual máximo 2.693e-8 relativo. Archivo durable pendiente. No formación prospectiva.

- E7.1a archivo CI 45595753850 SUCCESS, ZIP/4 NPZ/log/JSON verificados; blobs y tamaños releídos en 4ca90002. Gate instrumento aprobado, GO cerrar prerregistro E7.1b.

- E7.1b prerregistrada: 15 escenas, 4 pares prospectivos con puntos idénticos y orden distinto, 4 plantillas relajadas, 2 K0, replay; 7 épocas y criterios estructurales/funcionales fijos (>=3/4), límites de instrumento y negativos. Presupuesto15min gratuito. No datos de formación aún.

- E7.1b implementación y CI: 15 escenas, propietarios/recetas/NPZ por época, replay exacto, CDF por ventana, pares predefinidos; resultado negativo científico separado de fallo instrumental. Perfil heredado, suelo archivado como condición exógena, sin nueva ley. Presupuesto15min.

- E7.1b CI 45596252199: resultado NEGATIVO válido, 15/15 escenas y 105/105 ventanas, 0/4 pares cumplen, replay exacto, 2/4 diversidad estructural persistente sin función fuerte. Archivo pendiente. No repetir montículos. Usuario aprueba reorientación P2–P7 y representación material híbrida multiescala.

- E7.1b archivo 45604296564 SUCCESS, 150 outputs SHA/CRC y cotejo log/JSON verificados, blobs/tamaños releídos en 36d93f5d. Resultado negativo cerrado sin modificar gates. Dirección P2–P7 aprobada; ADR-WL-0001 y P2 consolidación completos documentalmente. GO P3.1: flujo real y contrato hidrológico, no repetir H11c cerrado ni confundir acumulación con descarga real.

- P3.1 prerregistro: inspección exacta OverlandFlow/Green-Ampt 2.11, q nativo m²/s pese a metadato, dx2, borde abierto, lluvia variable, stock subsuperficial finito, 5 escenarios y restart mismo host sin futuro. Sin erosión. Presupuesto5min.

- P3.1 código y workflow: flujo nativo de Almeida por enlace, balance independiente de frontera, Green-Ampt con capacidad finita, 5 casos y consumidor separado mismo host. Sin erosión ni acoplamiento nuevo. CI5min.

- P3.1 CI 45605745963 FAILURE instrumental antes de física, cinco iniciales roundtrip; máscara heredó subclase status_at_node. Corrección np.asarray preregistrada, sin cambios de gates. Archivo pendiente.

- Archivo primera incidencia P3.1 CI 45606846131 SUCCESS: 5 NC iniciales, ZIP/JSON/log cotejados, blobs/tamaños verificados. Corrección exclusiva máscara ndarray y retry mismo protocolo.

- P3.1 retry CI 45606950182 SUCCESS, 5×960 pasos, residual máximo1.4211e-14m³, export activo .702515m³, infiltración .355207m³, 480 sufijos y 5 controles exactos. Archivo pendiente. Candidato P3.2 RiverBedDynamics (no incluido en árbol v2.11.0, verificar distribución y contrato antes de adoptar).

- P3.1 archivo CI 45607198410 SUCCESS: 27 archivos internos cotejados, 20 outputs referencia, 2 consumidor, checkpoint y resúmenes; ZIP/log/JSON/CRC y blobs/tamaños verificados. P3.1 cerrada con límites. GO P3.2 evaluación proveedor conservativo.

- P3.2 selección: RBD candidato, MIT fork exacto, no reemplazar Landlab2.11. Sin límite de bedrock y outlet nativo anulado en fuente inspeccionada; comprobar con controles, no fingir export libre. P3.2a prerregistro: fixtures, live hidráulica→Exner unidireccional, stock finito y balances, 5min.

- P3.2a implementación y CI: recuperación hash exacto de 10 módulos y licencia, sin instalar fork; Exner nativo, controles flat/dry/export/depleted/replay, 3 propuestas desde hidrodinámica viva de terreno congelado, balances y stock finito; 5min.

- P3.2a CI 45615897598 SUCCESS: 9 escenarios (6 fixtures/replay y 3 live), gates vacíos, residual sólido max1.066e-15m³, depleted rechazado, export prescrito .00002m³, natural0. Live 120 pasos, RBD compatible. Archivo pendiente.

- P3.2a archivo 45619205788 SUCCESS; 12 outputs SHA/CRC y JSON/log verificados, blob/tamaño releído, gate aprobado limitado. P3.2b prerregistro de remapeo sólido/agua porosa con liberación local, capacidad dinámica y guard GA equivalente; no nueva física ni ajuste de gates.

- P3.2b implementación: 5 casos analíticos, 5 muestras Exner retenidas, 8 controles (incluye dominio GA y K0), operador conservativo con release proporcional a salidas brutas y capacidad derivada, GA nativo solo como tasa equivalente. CI5min; no acople vivo todavía.

- P3.2b CI45619612376 FAILURE: 5 analíticas, 8 controles, 4/5 muestras; error de comparación halo candidato bruto vs stock autoritativo (imposed_export). Corrección prospectiva de comparador core/owner + flujo halo explícito, sin modificar física/umbrales. Archivo pendiente.

- P3.2b archivo fallo 45624204131 SUCCESS, ZIP/9 outputs/JSON/log verificados; tamaños releídos. Retry con comparaciones explícitas de stock core/owner y flujo halo; sin cambiar operador físico ni tolerancias.

- P3.2b retry 45624348735 SUCCESS, 5 analíticas/5 muestras/8 controles, owner error5.55e-17m³, export/halo exacto, 9 NPZ previos idénticos. Archivo pendiente. GO P3.3 tras retención.

- P3.2b archivo 45625575751 SUCCESS, 10 outputs/JSON/log cotejados y 9 hashes históricos idénticos; tamaños verificados. P3.3 prerregistro: 240 pasos, 4 modos, ablación relieve→agua y bifurcación mismo estado paso121, balances y restart. Sin export natural sólido, sin nuevo solver. 5min.

- P3.3 implementación: 4 modos ×240 pasos, RBD y GA aislados por paso, hidráulica persistente, inventario único, bifurcación paso121, owner NC y consumidor mismo host, trazas completas de arrays. CI5min.

- P3.3 CI45625952812 FAILURE antes de física: HTTP429 descarga upstream. Sin artefacto. Archivar log y recuperar fuente exacta MIT como evidencia ZIP por blob API, no vendor activo; retry solo cambio de origen de bytes verificados.

- P3.3 fuente/fallo archivados CI45626242896 SUCCESS, ZIP38429bytes SHA256 b44242fead0857610603a75aec4148f3afa684f35307ca9005f9bc453a64720a; 11 miembros Git SHA1/tamaño y SHA256, log completo y blobs/tamaños verificados. Retry desde evidencia local exacta, no cambio de física.

- P3.3 retry CI45626406906 SUCCESS, 4×240 pasos, erosión/deposición acumuladas .0988295m³, Δz .0275477m, bifurcación q diff .00133745m²/s exacta rama física, restart120 y5 controles. Archivo pendiente, GO dictamen P3 y P4 tras retención.
