# Bitácora

- 2026-11-24: usuario aprobó diseño y migración. Origen inspeccionado, destino inicialmente vacío; corte 0448ea6e5e02d77142cec270ea99ac0773493fc1.
- Se recuperaron 12 commits desde dc1fdff. P0-D/P0-E completos, P1 aprobada, P1-A original negativa instrumental; P1-A-R CI 45488606315 y log de 12 casos con gates satisfechos. No se ejecutó física nueva. P1-B es siguiente fase científica.
- Se inició traspaso con workflow exclusivo de copia, comprobación de blobs y archivo de evidencia. Estado requiere lectura verificada y cobertura antes de avanzar.
- WL-M2: transferencia 45491750969 SUCCESS; 314/314 blobs exactos, 58/58 artefactos digest verificado. Origen permanece en 0448ea6e5e02d77142cec270ea99ac0773493fc1.
- WL-M3: archivo de logs 45492349703 SUCCESS, 481/487 ZIP, seis 404 sin jobs. Manifiestos y archivos recuperados en evidence/. Ningún artefacto inventariado faltante.
- WL-M4: verificación de árboles, evidencia, documentos, enlaces y simulación de continuidad (evidence/MIGRATION_VERIFICATION.md). Se archivan workflows de traspaso. P1-B.1 es siguiente unidad científica; STATE se actualizará al final.
- P1-B.1: inspección de contratos y fuentes Landlab v2.11.0; identificado déficit de flujos brutos en H33. Prerregistro de receptor por partición del suelo existente, entalpía incremental y transporte sedimentario. GO documental condicionado a B.2a; sin ejecución física nueva. Ver research/P1-B/P1-B.1_CONTRATO_Y_PRERREGISTRO.md.

- P1-B.2a: implementado adaptador de solo lectura de campos nativos SPACE (tasas × dt × área/densidad según unidades), comparador de hashes P0-D, controles y CI de 15 min. Referencias archivadas con digest verificado al ejecutar. Lanzamiento por push; resultado pendiente, no habilitar B.2b.
