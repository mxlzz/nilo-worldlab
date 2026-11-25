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
