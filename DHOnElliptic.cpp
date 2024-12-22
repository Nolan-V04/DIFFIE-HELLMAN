#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>

using namespace std;

struct BigInt {
    bool sign;  // true for positive, false for negative
    vector<uint32_t> BigIntDigits;
};

const uint32_t BASE = 1000000000; // Base for storing digits (9 decimal digits)

// Convert hex character to integer
uint32_t hexToInt(char hexChar) {
    if (hexChar >= '0' && hexChar <= '9') {
        return hexChar - '0';
    } else if (hexChar >= 'A' && hexChar <= 'F') {
        return hexChar - 'A' + 10;
    } else if (hexChar >= 'a' && hexChar <= 'f') {
        return hexChar - 'a' + 10;
    } else {
        throw invalid_argument("Invalid hex character");
    }
}

// Multiply BigInt by 16
void multiplyBy16(BigInt& bigInt) {
    uint64_t carry = 0;
    for (auto& digit : bigInt.BigIntDigits) {
        uint64_t product = static_cast<uint64_t>(digit) * 16 + carry;
        digit = product % BASE;
        carry = product / BASE;
    }
    if (carry > 0) {
        bigInt.BigIntDigits.push_back(static_cast<uint32_t>(carry));
    }
}

// Add a single digit to BigInt
void addDigit(BigInt& bigInt, uint32_t digit) {
    uint64_t carry = digit;
    for (auto& element : bigInt.BigIntDigits) {
        uint64_t sum = static_cast<uint64_t>(element) + carry;
        element = sum % BASE;
        carry = sum / BASE;
        if (carry == 0) return;
    }
    if (carry > 0) {
        bigInt.BigIntDigits.push_back(static_cast<uint32_t>(carry));
    }
}

// Convert hex string to BigInt
BigInt hexToBigInt(const string& hex) {
    BigInt result = {true, {0}};
    for (char hexChar : hex) {
        multiplyBy16(result);
        addDigit(result, hexToInt(hexChar));
    }
    return result;
}

// Compare two BigInts: -1 if a < b, 0 if a == b, 1 if a > b
int compareBigInts(const BigInt& a, const BigInt& b) {
    // Xét dấu trước
    if (a.sign != b.sign) return a.sign ? 1 : -1;

    // Nếu cùng dấu, so sánh phần giá trị
    int magnitudeComparison = 0;
    if (a.BigIntDigits.size() > b.BigIntDigits.size()) magnitudeComparison = 1;
    else if (a.BigIntDigits.size() < b.BigIntDigits.size()) magnitudeComparison = -1;
    else {
        for (size_t i = a.BigIntDigits.size(); i-- > 0;) {
            if (a.BigIntDigits[i] > b.BigIntDigits[i]) {
                magnitudeComparison = 1;
                break;
            }
            if (a.BigIntDigits[i] < b.BigIntDigits[i]) {
                magnitudeComparison = -1;
                break;
            }
        }
    }

    // Nếu là số âm, đảo kết quả so sánh
    return a.sign ? magnitudeComparison : -magnitudeComparison;
}

BigInt subtractBigInts(const BigInt& a, const BigInt& b);

BigInt addBigInts(const BigInt& a, const BigInt& b) {
    if (a.sign == b.sign) {
        // Nếu cùng dấu, thực hiện cộng
        BigInt result = {a.sign, {}};
        uint64_t carry = 0;
        size_t maxSize = max(a.BigIntDigits.size(), b.BigIntDigits.size());
        result.BigIntDigits.resize(maxSize);

        for (size_t i = 0; i < maxSize; ++i) {
            uint64_t sum = carry;
            if (i < a.BigIntDigits.size()) sum += a.BigIntDigits[i];
            if (i < b.BigIntDigits.size()) sum += b.BigIntDigits[i];
            result.BigIntDigits[i] = sum % BASE;
            carry = sum / BASE;
        }
        if (carry > 0) {
            result.BigIntDigits.push_back(static_cast<uint32_t>(carry));
        }
        return result;
    } else {
        // Nếu khác dấu, chuyển về phép trừ
        if (!a.sign) {
            // a âm, b dương => b - |a|
            return subtractBigInts(b, {true, a.BigIntDigits});
        } else {
            // a dương, b âm => a - |b|
            return subtractBigInts(a, {true, b.BigIntDigits});
        }
    }
}


BigInt subtractBigInts(const BigInt& a, const BigInt& b) {
    if (!a.sign && !b.sign) {
        // Nếu cả hai là số âm, đảo phép trừ: (-a) - (-b) => b - a
        return subtractBigInts(b, a);
    }

    if (a.sign != b.sign) {
        // Nếu khác dấu: a - (-b) => a + b
        return addBigInts(a, {a.sign, b.BigIntDigits});
    }

    // Nếu cùng dấu, kiểm tra giá trị lớn hơn
    if (compareBigInts(a, b) < 0) {
        BigInt result = subtractBigInts(b, a);
        result.sign = !a.sign; // Kết quả đảo dấu
        return result;
    }

    // Thực hiện trừ thông thường
    BigInt result = a;
    int64_t borrow = 0;

    for (size_t i = 0; i < b.BigIntDigits.size() || borrow > 0; ++i) {
        int64_t sub = borrow + (i < b.BigIntDigits.size() ? b.BigIntDigits[i] : 0);
        if (result.BigIntDigits[i] < sub) {
            result.BigIntDigits[i] += BASE - sub;
            borrow = 1;
        } else {
            result.BigIntDigits[i] -= sub;
            borrow = 0;
        }
    }

    // Xóa các chữ số 0 dư thừa
    while (result.BigIntDigits.size() > 1 && result.BigIntDigits.back() == 0) {
        result.BigIntDigits.pop_back();
    }

    return result;
}

BigInt multiplyBigInts(const BigInt& a, const BigInt& b) {
    if (a.BigIntDigits.empty() || b.BigIntDigits.empty()) return {true, {0}};
    BigInt result = {true, vector<uint32_t>(a.BigIntDigits.size() + b.BigIntDigits.size() + 1, 0)};

    result.sign = (a.sign == b.sign);

    for (size_t i = 0; i < a.BigIntDigits.size(); ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < b.BigIntDigits.size(); ++j) {
            uint64_t product = static_cast<uint64_t>(a.BigIntDigits[i]) * b.BigIntDigits[j] + result.BigIntDigits[i + j] + carry;
            result.BigIntDigits[i + j] = product % BASE;
            carry = product / BASE;
        }
        if (carry > 0) {
            result.BigIntDigits[i + b.BigIntDigits.size()] += static_cast<uint32_t>(carry);
        }
    }

    while (result.BigIntDigits.size() > 1 && result.BigIntDigits.back() == 0) {
        result.BigIntDigits.pop_back();
    }

    return result;
}


BigInt intToBigInt(uint64_t value) {
    BigInt result = {true, {}};
    while (value > 0) {
        uint32_t digit = value % BASE;
        result.BigIntDigits.push_back(digit);
        value /= BASE;
    }
    if (result.BigIntDigits.empty()) {
        result.BigIntDigits.push_back(0);
    }
    
    return result;
}

pair<BigInt, BigInt> divideBigInts(const BigInt& dividend, const BigInt& divisor) {
    if (divisor.BigIntDigits.empty() || (divisor.BigIntDigits.size() == 1 && divisor.BigIntDigits[0] == 0)) {
        throw invalid_argument("Division by zero");
    }

    if (compareBigInts(dividend, divisor) < 0) {
        return {{true, {0}}, dividend}; // Thương = 0, Dư = dividend
    }

    BigInt quotient = {dividend.sign == divisor.sign, vector<uint32_t>(dividend.BigIntDigits.size(), 0)};
    BigInt remainder = {dividend.sign, dividend.BigIntDigits};

    size_t shift = dividend.BigIntDigits.size() - divisor.BigIntDigits.size();
    BigInt shiftedDivisor = divisor;
    shiftedDivisor.BigIntDigits.insert(shiftedDivisor.BigIntDigits.begin(), shift, 0);

    for (size_t i = shift + 1; i-- > 0;) {
        uint32_t q = 0;
        uint32_t low = 0, high = BASE - 1;

        while (low <= high) {
            uint32_t mid = low + (high - low) / 2;
            BigInt testProduct = multiplyBigInts(shiftedDivisor, {true, {mid}});

            if (compareBigInts(testProduct, remainder) <= 0) {
                q = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }

        quotient.BigIntDigits[i] = q;
        remainder = subtractBigInts(remainder, multiplyBigInts(shiftedDivisor, {true, {q}}));

        if (i > 0) {
            shiftedDivisor.BigIntDigits.erase(shiftedDivisor.BigIntDigits.begin());
        }
    }

    // Xóa các chữ số 0 dư thừa
    while (quotient.BigIntDigits.size() > 1 && quotient.BigIntDigits.back() == 0) {
        quotient.BigIntDigits.pop_back();
    }
    while (remainder.BigIntDigits.size() > 1 && remainder.BigIntDigits.back() == 0) {
        remainder.BigIntDigits.pop_back();
    }

    return {quotient, remainder};
}

void printBigInt(const BigInt& bigInt) {
    if (bigInt.BigIntDigits.empty()) {
        cout << "0\n";
        return;
    }
    if (!bigInt.sign) cout << "-";

    // Print the most significant digit without leading zeros
    cout << bigInt.BigIntDigits.back();

    // Print the remaining digits with leading zeros
    for (int i = bigInt.BigIntDigits.size() - 2; i >= 0; --i) {
        cout << setw(9) << setfill('0') << bigInt.BigIntDigits[i];
    }
    cout << endl;
}

// Hàm modBigInts
BigInt modBigInts(const BigInt& a, const BigInt& b) {
    if (b.BigIntDigits.empty() || (b.BigIntDigits.size() == 1 && b.BigIntDigits[0] == 0)) {
        throw invalid_argument("Modulo by zero");
    }

    // Calculate the remainder
    BigInt remainder = divideBigInts(a, b).second;

    // If the remainder is negative, adjust it by adding b to it
    while (!remainder.sign) {
        remainder = addBigInts(remainder, b);
    }

    // Ensure the remainder is non-negative and smaller than b
    if (compareBigInts(remainder, b) >= 0) {
        remainder = subtractBigInts(remainder, b);
    }

    return remainder;
}


BigInt modExp(const BigInt& base, const BigInt& exponent, const BigInt& modulus) {
    if (modulus.BigIntDigits.empty() || (modulus.BigIntDigits.size() == 1 && modulus.BigIntDigits[0] == 0)) {
        throw invalid_argument("Modulo by zero");
    }

    BigInt result = {true, {1}}; // Initialize result as 1
    BigInt baseMod = modBigInts(base, modulus); // Ensure base is within modulus range
    BigInt exp = exponent; // Copy of the exponent for modification

    while (!exp.BigIntDigits.empty() && !(exp.BigIntDigits.size() == 1 && exp.BigIntDigits[0] == 0)) {
        // Check if the current exponent is odd
        if (exp.BigIntDigits[0] % 2 == 1) {
            result = modBigInts(multiplyBigInts(result, baseMod), modulus);
        }

        // Divide the exponent by 2
        exp = divideBigInts(exp, {true, {2}}).first;

        // Square the base and take modulo
        baseMod = modBigInts(multiplyBigInts(baseMod, baseMod), modulus);
    }

    return result;
}

// Extended Euclidean Algorithm to find gcd and coefficients for BigInt
tuple<BigInt, BigInt, BigInt> extendedGCD(const BigInt& a, const BigInt& b) {
    // Base case
    if (compareBigInts(b, {true, {0}}) == 0) {
        return make_tuple(a, BigInt{true, {1}}, BigInt{true, {0}});
    }

tuple<BigInt, BigInt, BigInt> result = extendedGCD(b, modBigInts(a, b));  // Recursive call
    BigInt gcd = get<0>(result);
    BigInt x1 = get<1>(result);
    BigInt y1 = get<2>(result);

    // Update coefficients
    BigInt x = y1;
    BigInt quotient = divideBigInts(a, b).first;
    BigInt y = subtractBigInts(x1, multiplyBigInts(quotient, y1));  // Assume BigInt division and subtraction are implemented
    
    return make_tuple(gcd, x, y);
}

// Modular Inverse for BigInt
BigInt modInverse(const BigInt& a, const BigInt& p) {
    BigInt t = {true, {0}};  // t = 0
    BigInt newT = {true, {1}};  // newT = 1
    BigInt r = p;  // r = p
    BigInt newR = a;  // newR = a

    while (compareBigInts(newR, {true, {0}}) != 0) {  // while newR != 0
        BigInt quotient = divideBigInts(r, newR).first;  // quotient = r / newR

        // Update t and newT
        BigInt tempT = t;
        t = newT;
        newT = subtractBigInts(tempT, multiplyBigInts(quotient, newT));

        // Update r and newR
        BigInt tempR = r;
        r = newR;
        newR = subtractBigInts(tempR, multiplyBigInts(quotient, newR));
    }

    // Ensure gcd(a, p) == 1
    if (compareBigInts(r, {true, {1}}) != 0) {
        throw invalid_argument("a and p are not coprime, so no modular inverse exists.");
    }

    // Adjust t to be positive
    if (!t.sign) {
        t = addBigInts(t, p);
    }

    return modBigInts(t, p);  // Ensure result is in the range [0, p - 1]
}


// Add two points on elliptic
pair<BigInt, BigInt> addPoints(const BigInt& x1, const BigInt& y1, const BigInt& x2, const BigInt& y2, const BigInt& a, const BigInt& p) {
    if (compareBigInts(x1, x2) == 0 && compareBigInts(y1, y2) == 0) { // P = Q
        //lambda = (3 * x1^2 + a) / (2 * y1) mod p
        BigInt numerator = addBigInts(multiplyBigInts(intToBigInt(3), multiplyBigInts(x1, x1)), a);
        BigInt denominator = multiplyBigInts(intToBigInt(2), y1);
        BigInt lambda = modBigInts(multiplyBigInts(numerator, modInverse(denominator, p)), p);
        
        //xR = lambda^2 - 2 * x1 mod p
        BigInt xR = modBigInts(subtractBigInts(multiplyBigInts(lambda, lambda), multiplyBigInts(intToBigInt(2), x1)), p);
       
        //yR = lambda * (x1 - xR) - y1 mod p
        BigInt yR = modBigInts(subtractBigInts(multiplyBigInts(lambda, subtractBigInts(x1, xR)), y1), p);
        
        return {xR, yR};
    } else { // P != Q
        //lambda = (y2 - y1) / (x2 - x1) mod p
        BigInt numerator = subtractBigInts(y2, y1);
        BigInt denominator = subtractBigInts(x2, x1);
        BigInt lambda = modBigInts(multiplyBigInts(numerator, modInverse(denominator, p)), p);

        //xR = lambda^2 - x1 - x2 mod p
        BigInt xR = modBigInts(subtractBigInts(subtractBigInts(multiplyBigInts(lambda, lambda), x1), x2), p);

        //yR = lambda * (x1 - xR) - y1 mod p
        BigInt yR = modBigInts(subtractBigInts(multiplyBigInts(lambda, subtractBigInts(x1, xR)), y1), p);

        return {xR, yR};
    }

}

string BigIntToHex(const BigInt& bigInt) {
    if (bigInt.BigIntDigits.empty()) {
        return "0"; // Handle empty BigInt as zero
    }

    string hexString;

    // Traverse digits in reverse (to handle big-endian output)
    for (int i = bigInt.BigIntDigits.size() - 1; i >= 0; --i) {
        stringstream ss;
        ss << hex << uppercase << bigInt.BigIntDigits[i];

        string hexPart = ss.str();

        // Ensure each part is padded to 9 hex digits (except the most significant one)
        if (i != bigInt.BigIntDigits.size() - 1) {
            while (hexPart.length() < 9) {
                hexPart = "0" + hexPart;
            }
        }

        hexString += hexPart;
    }

    // Remove leading zeros in the final hexadecimal string
    size_t firstNonZero = hexString.find_first_not_of('0');
    if (firstNonZero != string::npos) {
        hexString = hexString.substr(firstNonZero);
    } else {
        hexString = "0"; // Handle case where the number is zero
    }

    return hexString;
}

void readInput(const string& filename, BigInt& p, BigInt& a, BigInt& b, BigInt& xP, BigInt& yP, BigInt& xQ, BigInt& yQ) {
    ifstream fin(filename);
    if(!fin.is_open()) {
        cout << "Cannot open file" + filename + "! Retry." << endl;
    }
    string pHex, aHex, bHex, xPHex, yPHex, xQHex, yQHex;
    fin >> pHex >> aHex >> bHex >> xPHex >> yPHex >> xQHex >> yQHex;

    p = hexToBigInt(pHex);
    a = hexToBigInt(aHex);
    b = hexToBigInt(bHex);
    xP = hexToBigInt(xPHex);
    yP = hexToBigInt(yPHex);
    xQ = hexToBigInt(xQHex);
    yQ = hexToBigInt(yQHex);

    fin.close();
}

void writeOutput(const string& filename, BigInt xR, BigInt yR) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cout << "Cannot open file" + filename + "! Retry." << endl;
        exit(1);
    }
    string xRHex = BigIntToHex(xR);
    string yRHex = BigIntToHex(yR);
    fout << xRHex << " " << yRHex << endl;
    fout.close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    auto start = chrono::high_resolution_clock::now();

    BigInt p, a, b, xP, yP, xQ, yQ;
    
    try {
        // Đọc dữ liệu từ file
        readInput(argv[1], p, a, b, xP, yP, xQ, yQ);

        // Tính toán điểm cộng
        pair<BigInt, BigInt> result = addPoints(xP, yP, xQ, yQ, a, p);
        BigInt xR = result.first;
        BigInt yR = result.second;

        // Ghi kết quả ra file
        writeOutput(argv[2], xR, yR);

    } catch (const exception& e) {
        cout << e.what() << endl;
        return 1;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Execution time: " << elapsed.count() << " seconds" << endl;


    return 0;
}