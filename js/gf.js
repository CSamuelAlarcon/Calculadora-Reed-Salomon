/**
 * Aritmética de cuerpos de Galois GF(2^m) para la calculadora Reed-Solomon.
 *
 * Cada símbolo del código es un elemento de GF(2^m), representado como un
 * entero de 0 a 2^m - 1. Suma y resta son XOR; multiplicación y división se
 * resuelven con tablas de logaritmos/antilogaritmos construidas a partir de
 * un polinomio primitivo (ver TFI, sección 4.1 del anteproyecto).
 */

// Polinomios primitivos estándar por tamaño de símbolo m (en binario, con el
// bit líder x^m implícito incluido para poder hacer XOR directo al reducir).
const PRIMITIVE_POLYNOMIALS = {
  3: 0b1011, // x^3 + x + 1
  4: 0b10011, // x^4 + x + 1        (usado en el ejemplo del anteproyecto, GF(16))
  5: 0b100101, // x^5 + x^2 + 1
  6: 0b1000011, // x^6 + x + 1
  7: 0b10001001, // x^7 + x^3 + 1
  8: 0b100011101, // x^8 + x^4 + x^3 + x^2 + 1  (estándar RS(255,k), QR, AES)
};

class GaloisField {
  /**
   * @param {number} m - tamaño de símbolo en bits (3 a 8 recomendado)
   * @param {number} [primitivePoly] - polinomio primitivo en binario;
   *   si se omite, se usa el estándar de PRIMITIVE_POLYNOMIALS[m]
   */
  constructor(m, primitivePoly = PRIMITIVE_POLYNOMIALS[m]) {
    if (!Number.isInteger(m) || m < 2 || m > 16) {
      throw new Error(`m debe ser un entero entre 2 y 16 (recibido: ${m})`);
    }
    if (!primitivePoly) {
      throw new Error(
        `No hay polinomio primitivo por defecto para m=${m}. Debe indicarse manualmente.`
      );
    }

    this.m = m;
    this.size = 1 << m; // 2^m: cantidad total de elementos (incluye el 0)
    this.max = this.size - 1; // 2^m - 1: elementos no nulos (orden del grupo multiplicativo)
    this.primitivePoly = primitivePoly;

    // expTable se duplica (largo 2*max) para poder sumar exponentes en la
    // multiplicación sin aplicar % max en cada operación.
    this.expTable = new Array(this.max * 2);
    this.logTable = new Array(this.size);

    this._buildTables();
  }

  _buildTables() {
    let x = 1;
    for (let i = 0; i < this.max; i++) {
      this.expTable[i] = x;
      this.logTable[x] = i;

      x <<= 1; // multiplicar por alpha = desplazar 1 bit a la izquierda
      if (x & this.size) {
        // se "desbordó" el bit m: reducir aplicando la regla del polinomio primitivo
        x ^= this.primitivePoly;
      }
    }
    for (let i = this.max; i < this.max * 2; i++) {
      this.expTable[i] = this.expTable[i - this.max];
    }
    this.logTable[0] = undefined; // log(0) no está definido
  }

  /** Suma en GF(2^m): equivale a XOR. */
  add(a, b) {
    return a ^ b;
  }

  /** Resta en GF(2^m): idéntica a la suma (es XOR). */
  sub(a, b) {
    return a ^ b;
  }

  /** Multiplicación vía tabla de logaritmos: log(a)+log(b) -> antilog. */
  mul(a, b) {
    this._checkElement(a);
    this._checkElement(b);
    if (a === 0 || b === 0) return 0;
    return this.expTable[this.logTable[a] + this.logTable[b]];
  }

  /** División vía tabla de logaritmos: log(a)-log(b) -> antilog. */
  div(a, b) {
    this._checkElement(a);
    this._checkElement(b);
    if (b === 0) throw new Error('División por cero en GF(2^m)');
    if (a === 0) return 0;
    return this.expTable[(this.logTable[a] - this.logTable[b] + this.max) % this.max];
  }

  /** Potencia a^n dentro del cuerpo (n puede ser negativo). */
  pow(a, n) {
    this._checkElement(a);
    if (a === 0) return n === 0 ? 1 : 0;
    let e = (this.logTable[a] * n) % this.max;
    if (e < 0) e += this.max;
    return this.expTable[e];
  }

  /** Inverso multiplicativo de a (a != 0). */
  inv(a) {
    this._checkElement(a);
    if (a === 0) throw new Error('El elemento 0 no tiene inverso multiplicativo');
    return this.expTable[this.max - this.logTable[a]];
  }

  /** alpha^exp, aceptando exponentes negativos o mayores a max (se reducen mod max). */
  alphaTo(exp) {
    let e = exp % this.max;
    if (e < 0) e += this.max;
    return this.expTable[e];
  }

  /** log_alpha(a). */
  log(a) {
    this._checkElement(a);
    if (a === 0) throw new Error('log(0) no está definido en GF(2^m)');
    return this.logTable[a];
  }

  _checkElement(a) {
    if (!Number.isInteger(a) || a < 0 || a > this.max) {
      throw new Error(`Elemento fuera de rango para GF(2^${this.m}): ${a}`);
    }
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { GaloisField, PRIMITIVE_POLYNOMIALS };
}
