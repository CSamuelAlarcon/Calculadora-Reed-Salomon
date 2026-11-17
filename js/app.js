/**
 * Interfaz de la calculadora Reed-Solomon: conecta los módulos gf.js,
 * rsEncoder.js y rsDecoder.js con los formularios de la página.
 */
(function () {
  'use strict';

  const el = (id) => document.getElementById(id);

  // --- Formato de exponentes como superíndice real (x³, GF(2⁸), etc.) ---
  const SUPERSCRIPT_DIGITS = {
    '0': '⁰',
    '1': '¹',
    '2': '²',
    '3': '³',
    '4': '⁴',
    '5': '⁵',
    '6': '⁶',
    '7': '⁷',
    '8': '⁸',
    '9': '⁹',
  };
  function toSuperscript(n) {
    return String(n)
      .split('')
      .map((ch) => SUPERSCRIPT_DIGITS[ch] || ch)
      .join('');
  }

  // --- Referencias del DOM ---
  const mSelect = el('mSelect');
  const kInput = el('kInput');
  const tInput = el('tInput');
  const configSummary = el('configSummary');
  const configError = el('configError');
  const loadExampleBtn = el('loadExampleBtn');

  const messageInput = el('messageInput');
  const encodeBtn = el('encodeBtn');
  const encodeError = el('encodeError');
  const encodeResults = el('encodeResults');
  const logTableContainer = el('logTableContainer');
  const generatorOutput = el('generatorOutput');
  const codewordOutput = el('codewordOutput');

  const channelSection = el('channel-section');
  const tDisplay = el('tDisplay');
  const receivedSymbolsEl = el('receivedSymbols');
  const decodeBtn = el('decodeBtn');

  const decodeSection = el('decode-section');
  const syndromesOutput = el('syndromesOutput');
  const locatorOutput = el('locatorOutput');
  const errorsTableBody = document.querySelector('#errorsTable tbody');
  const correctedOutput = el('correctedOutput');
  const decodeStatus = el('decodeStatus');

  // --- Estado de la calculadora ---
  let gf = null;
  let currentM = null;
  let currentNsym = null;
  let receivedSymbols = null;

  // --- Configuración inicial: poblar el selector de m ---
  for (let m = 3; m <= 8; m++) {
    const opt = document.createElement('option');
    opt.value = String(m);
    opt.textContent = 'm = ' + m + ' (GF(' + (1 << m) + '))';
    if (m === 4) opt.selected = true;
    mSelect.appendChild(opt);
  }

  function currentConfig() {
    const m = parseInt(mSelect.value, 10);
    const k = parseInt(kInput.value, 10);
    const t = parseInt(tInput.value, 10);
    const nsym = 2 * t;
    const n = k + nsym;
    const maxN = (1 << m) - 1;
    return { m, k, t, nsym, n, maxN };
  }

  function renderConfigSummary() {
    const { m, k, t, nsym, n, maxN } = currentConfig();

    if (!Number.isInteger(k) || k < 1 || !Number.isInteger(t) || t < 1) {
      configSummary.textContent = '';
      configError.textContent = 'k y t deben ser enteros mayores o iguales a 1.';
      return;
    }
    if (n > maxN) {
      configSummary.textContent = '';
      configError.textContent =
        'n = k + 2t = ' + n + ' supera el máximo permitido para GF(2' + toSuperscript(m) + '), que es ' + maxN +
        '. Reducí k o t, o aumentá m.';
      return;
    }
    configError.textContent = '';
    configSummary.textContent =
      'Código RS(' + n + ', ' + k + ') sobre GF(' + (1 << m) + ') — corrige hasta t = ' + t +
      ' símbolos erróneos por bloque (2t = ' + nsym + ' símbolos de paridad).';
  }

  [mSelect, kInput, tInput].forEach((node) => node.addEventListener('input', renderConfigSummary));
  renderConfigSummary();

  loadExampleBtn.addEventListener('click', () => {
    mSelect.value = '4';
    kInput.value = '3';
    tInput.value = '2';
    messageInput.value = '12,6,3';
    renderConfigSummary();
    encodeMessage();
  });

  /** Convierte "12,6,3" en [12,6,3], validando cada símbolo contra el máximo del cuerpo. */
  function parseSymbolList(text, max) {
    const parts = text
      .split(',')
      .map((s) => s.trim())
      .filter((s) => s.length > 0);
    if (parts.length === 0) {
      throw new Error('Ingresá al menos un símbolo.');
    }
    const symbols = [];
    for (const part of parts) {
      const value = Number(part);
      if (!Number.isInteger(value) || value < 0 || value > max) {
        throw new Error('"' + part + '" no es válido: debe ser un entero entre 0 y ' + max + '.');
      }
      symbols.push(value);
    }
    return symbols;
  }

  /** Formatea un polinomio LEADING-FIRST (índice 0 = mayor grado) como texto legible. */
  function formatPolynomial(coeffsLeadingFirst) {
    const degree = coeffsLeadingFirst.length - 1;
    const terms = [];
    coeffsLeadingFirst.forEach((c, idx) => {
      if (c === 0) return;
      const exp = degree - idx;
      const coefPrefix = c === 1 && exp > 0 ? '' : String(c); // omitir el "1" delante de x, x², etc.
      if (exp === 0) terms.push(String(c));
      else if (exp === 1) terms.push(coefPrefix + 'x');
      else terms.push(coefPrefix + 'x' + toSuperscript(exp));
    });
    return terms.length ? terms.join(' + ') : '0';
  }

  function renderLogTable(gfInstance) {
    let html = '<table><thead><tr><th>i</th><th>&alpha;<sup>i</sup></th></tr></thead><tbody>';
    for (let i = 0; i < gfInstance.max; i++) {
      html += '<tr><td>' + i + '</td><td>' + gfInstance.alphaTo(i) + '</td></tr>';
    }
    html += '</tbody></table>';
    logTableContainer.innerHTML = html;
  }

  /**
   * Dibuja una fila de símbolos (codeword). Si editable=true, cada símbolo
   * se muestra como <input> y se llama a onChange(indice, inputElement)
   * cuando el usuario lo modifica.
   */
  function renderSymbolRow(container, symbols, editable, onChange) {
    container.innerHTML = '';
    const n = symbols.length;
    symbols.forEach((value, idx) => {
      const wrapper = document.createElement('div');
      wrapper.className = 'symbol';

      const label = document.createElement('span');
      label.className = 'symbol-exp';
      const exponent = n - 1 - idx; // leading-first: índice 0 = x^(n-1)
      label.textContent = 'x' + toSuperscript(exponent);
      wrapper.appendChild(label);

      if (editable) {
        const input = document.createElement('input');
        input.type = 'number';
        input.min = '0';
        input.value = String(value);
        input.dataset.index = String(idx);
        input.addEventListener('change', () => onChange(idx, input));
        wrapper.appendChild(input);
      } else {
        const span = document.createElement('span');
        span.className = 'symbol-value';
        span.textContent = String(value);
        wrapper.appendChild(span);
      }

      container.appendChild(wrapper);
    });
  }

  function encodeMessage() {
    encodeError.textContent = '';
    decodeSection.classList.add('hidden');
    channelSection.classList.add('hidden');
    encodeResults.classList.add('hidden');

    const { m, k, nsym, n, maxN, t } = currentConfig();
    if (n > maxN || k < 1 || t < 1) {
      encodeError.textContent = 'Revisá la configuración del código antes de codificar (sección 1).';
      return;
    }

    let symbols;
    try {
      symbols = parseSymbolList(messageInput.value, (1 << m) - 1);
    } catch (e) {
      encodeError.textContent = e.message;
      return;
    }
    if (symbols.length !== k) {
      encodeError.textContent =
        'Ingresaste ' + symbols.length + ' símbolo(s), pero k = ' + k + '. Ajustá el mensaje o el valor de k.';
      return;
    }

    gf = new GaloisField(m);
    currentM = m;
    currentNsym = nsym;

    const { codeword, generator } = rsEncode(gf, symbols, nsym);
    receivedSymbols = codeword.slice();

    renderLogTable(gf);
    generatorOutput.textContent = 'G(x) = ' + formatPolynomial(generator);
    renderSymbolRow(codewordOutput, codeword, false);
    encodeResults.classList.remove('hidden');

    tDisplay.textContent = String(t);
    renderSymbolRow(receivedSymbolsEl, receivedSymbols, true, (idx, input) => {
      const val = Number(input.value);
      const max = (1 << currentM) - 1;
      if (!Number.isInteger(val) || val < 0 || val > max) {
        input.value = String(receivedSymbols[idx]);
        return;
      }
      receivedSymbols[idx] = val;
    });
    channelSection.classList.remove('hidden');
  }

  encodeBtn.addEventListener('click', encodeMessage);

  decodeBtn.addEventListener('click', () => {
    if (!gf || !receivedSymbols) return;

    const result = rsDecode(gf, receivedSymbols, currentNsym);

    syndromesOutput.textContent = '[' + result.syndromes.join(', ') + ']';

    // locator viene en orden ascendente (Λ[0] = término independiente);
    // se invierte para reutilizar formatPolynomial, que espera leading-first.
    locatorOutput.textContent = 'Λ(x) = ' + formatPolynomial(result.locator.slice().reverse());

    errorsTableBody.innerHTML = '';
    result.errors.forEach((err) => {
      const tr = document.createElement('tr');
      const tdPos = document.createElement('td');
      tdPos.textContent = 'x' + toSuperscript(err.exponent);
      const tdMag = document.createElement('td');
      tdMag.textContent = String(err.magnitude);
      tr.appendChild(tdPos);
      tr.appendChild(tdMag);
      errorsTableBody.appendChild(tr);
    });

    renderSymbolRow(correctedOutput, result.corrected, false);

    decodeStatus.classList.remove('status-ok', 'status-fail');
    if (!result.hadErrors) {
      decodeStatus.classList.add('status-ok');
      decodeStatus.textContent = 'No se detectaron errores: el codeword recibido es válido.';
    } else if (result.success) {
      decodeStatus.classList.add('status-ok');
      decodeStatus.textContent =
        'Se detectaron y corrigieron ' + result.numErrorsFound + ' error(es) correctamente.';
    } else {
      decodeStatus.classList.add('status-fail');
      decodeStatus.textContent = result.reason || 'No se pudo corregir el codeword recibido.';
    }

    decodeSection.classList.remove('hidden');
  });
})();
