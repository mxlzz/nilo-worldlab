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
