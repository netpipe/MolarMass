#include <QApplication>
#include <QWidget>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QVBoxLayout>
#include <QHash>
#include <QString>
#include <QStringList>
#include <QChar>
#include <QObject>

class FormulaParser {
public:
    static bool parse(const QString &input, double &mass, QString &error)
    {
        QString f = input.simplified();
        f.remove(QChar(' '));

        // Different dot-like separators used in hydrates/adducts.
        f.replace(QChar(0x00B7), QChar('.')); // middle dot ·
        f.replace(QChar(0x22C5), QChar('.')); // dot operator ⋅
        f.replace(QChar(0x2022), QChar('.')); // bullet •
        f.replace(QChar('*'), QChar('.'));
        f.replace(QChar(':'), QChar('.'));

        if (f.isEmpty()) {
            error = QStringLiteral("Empty formula.");
            return false;
        }

        QStringList parts = f.split(QChar('.'), QString::SkipEmptyParts);
        if (parts.isEmpty()) {
            error = QStringLiteral("Empty formula.");
            return false;
        }

        double total = 0.0;

        for (int i = 0; i < parts.size(); ++i) {
            QString seg = parts.at(i).trimmed();
            if (seg.isEmpty())
                continue;

            int pos = 0;
            qint64 coeff = 1;

            // Optional leading coefficient.
            // For the first segment, leading coefficients are ignored because
            // molar mass is usually per formula unit. For hydrate segments after
            // a dot, coefficients are used: CuSO4.5H2O.
            if (pos < seg.size() && seg.at(pos).isDigit()) {
                int start = pos;
                while (pos < seg.size() && seg.at(pos).isDigit())
                    ++pos;

                bool validStart =
                    pos < seg.size() &&
                    (seg.at(pos).isUpper() || isOpening(seg.at(pos)));

                if (validStart) {
                    if (i > 0) {
                        bool ok = false;
                        long long v = seg.mid(start, pos - start).toLongLong(&ok);
                        if (!ok || v < 0) {
                            error = QStringLiteral("Invalid coefficient.");
                            return false;
                        }
                        coeff = static_cast<qint64>(v);
                    } else {
                        // Ignore leading coefficient on the first formula segment.
                    }
                } else {
                    pos = 0;
                }
            }

            if (pos >= seg.size()) {
                error = QStringLiteral("Missing formula.");
                return false;
            }

            double segmentMass = 0.0;
            if (!parseFormula(seg, pos, QChar(), segmentMass, error))
                return false;

            if (pos != seg.size()) {
                error = QStringLiteral("Unexpected character at position %1.").arg(pos);
                return false;
            }

            total += static_cast<double>(coeff) * segmentMass;
        }

        if (total <= 0.0) {
            error = QStringLiteral("Formula has no calculable mass.");
            return false;
        }

        mass = total;
        return true;
    }

private:
    static bool isOpening(QChar c)
    {
        return c == QChar('(') || c == QChar('[') || c == QChar('{');
    }

    static QChar matchingClosing(QChar c)
    {
        if (c == QChar('('))
            return QChar(')');
        if (c == QChar('['))
            return QChar(']');
        return QChar('}');
    }

    static bool parseFormula(const QString &s,
                             int &pos,
                             QChar closing,
                             double &mass,
                             QString &error)
    {
        mass = 0.0;
        const QHash<QString, double> &table = atomicMasses();

        while (pos < s.length()) {
            QChar c = s.at(pos);

            if (!closing.isNull() && c == closing)
                break;

            if (isOpening(c)) {
                QChar close = matchingClosing(c);
                ++pos;

                double inner = 0.0;
                if (!parseFormula(s, pos, close, inner, error))
                    return false;

                if (pos >= s.length() || s.at(pos) != close) {
                    error = QStringLiteral("Missing '%1' at position %2.")
                                .arg(QString(close))
                                .arg(pos);
                    return false;
                }

                ++pos;

                qint64 mult = 1;
                if (!parseMultiplier(s, pos, mult, error, true))
                    return false;

                if (inner <= 0.0) {
                    error = QStringLiteral("Empty group.");
                    return false;
                }

                mass += static_cast<double>(mult) * inner;

            } else if (c.isUpper()) {
                QString symbol;
                symbol.append(c);
                ++pos;

                while (pos < s.length() && s.at(pos).isLower()) {
                    symbol.append(s.at(pos));
                    ++pos;
                }

                qint64 mult = 1;
                if (!parseMultiplier(s, pos, mult, error, false))
                    return false;

                QHash<QString, double>::const_iterator it = table.constFind(symbol);
                if (it == table.constEnd()) {
                    error = QStringLiteral("Unknown element: %1.").arg(symbol);
                    return false;
                }

                mass += it.value() * static_cast<double>(mult);

            } else if (c == QChar('^')) {
                // Skip charge notation like ^2+, ^-, ^+2.
                ++pos;
                while (pos < s.length() &&
                       (s.at(pos).isDigit() ||
                        s.at(pos) == QChar('+') ||
                        s.at(pos) == QChar('-'))) {
                    ++pos;
                }

            } else if (c == QChar('+') || c == QChar('-')) {
                // Skip simple charge sign, possibly followed by digits: +2, -3.
                ++pos;
                while (pos < s.length() && s.at(pos).isDigit())
                    ++pos;

            } else if (c.isSpace()) {
                ++pos;

            } else if (c.isDigit()) {
                error = QStringLiteral("Unexpected number at position %1.").arg(pos);
                return false;

            } else {
                error = QStringLiteral("Unexpected character '%1' at position %2.")
                            .arg(QString(c))
                            .arg(pos);
                return false;
            }
        }

        return true;
    }

    static bool parseMultiplier(const QString &s,
                                int &pos,
                                qint64 &mult,
                                QString &error,
                                bool afterGroup)
    {
        mult = 1;

        if (pos >= s.length() || !s.at(pos).isDigit())
            return true;

        int start = pos;
        while (pos < s.length() && s.at(pos).isDigit())
            ++pos;

        int after = pos;

        // For groups, a number followed by + or - is usually a charge:
        // [Fe(CN)6]4-
        if (afterGroup &&
            after < s.length() &&
            (s.at(after) == QChar('+') || s.at(after) == QChar('-'))) {
            pos = after + 1;

            // Also skip any digits after the sign, for unusual charge forms.
            while (pos < s.length() && s.at(pos).isDigit())
                ++pos;

            mult = 1;
            return true;
        }

        bool ok = false;
        long long v = s.mid(start, pos - start).toLongLong(&ok);

        if (!ok || v < 0 || v > 1000000000LL) {
            error = QStringLiteral("Invalid multiplier at position %1.").arg(start);
            return false;
        }

        mult = static_cast<qint64>(v);
        return true;
    }

    static const QHash<QString, double> &atomicMasses()
    {
        // Standard atomic weights where available.
        // Radioactive elements use approximate masses.
        // D and T are included as common isotope symbols.
        static const QHash<QString, double> masses = {
            { QStringLiteral("H"),  1.008 },
            { QStringLiteral("D"),  2.01410177812 },
            { QStringLiteral("T"),  3.0160492779 },
            { QStringLiteral("He"), 4.002602 },
            { QStringLiteral("Li"), 6.94 },
            { QStringLiteral("Be"), 9.0121831 },
            { QStringLiteral("B"),  10.81 },
            { QStringLiteral("C"),  12.011 },
            { QStringLiteral("N"),  14.007 },
            { QStringLiteral("O"),  15.999 },
            { QStringLiteral("F"),  18.998403163 },
            { QStringLiteral("Ne"), 20.1797 },
            { QStringLiteral("Na"), 22.98976928 },
            { QStringLiteral("Mg"), 24.305 },
            { QStringLiteral("Al"), 26.9815385 },
            { QStringLiteral("Si"), 28.085 },
            { QStringLiteral("P"),  30.973761998 },
            { QStringLiteral("S"),  32.06 },
            { QStringLiteral("Cl"), 35.45 },
            { QStringLiteral("Ar"), 39.948 },
            { QStringLiteral("K"),  39.0983 },
            { QStringLiteral("Ca"), 40.078 },
            { QStringLiteral("Sc"), 44.955908 },
            { QStringLiteral("Ti"), 47.867 },
            { QStringLiteral("V"),  50.9415 },
            { QStringLiteral("Cr"), 51.9961 },
            { QStringLiteral("Mn"), 54.938044 },
            { QStringLiteral("Fe"), 55.845 },
            { QStringLiteral("Co"), 58.933194 },
            { QStringLiteral("Ni"), 58.6934 },
            { QStringLiteral("Cu"), 63.546 },
            { QStringLiteral("Zn"), 65.38 },
            { QStringLiteral("Ga"), 69.723 },
            { QStringLiteral("Ge"), 72.630 },
            { QStringLiteral("As"), 74.921595 },
            { QStringLiteral("Se"), 78.971 },
            { QStringLiteral("Br"), 79.904 },
            { QStringLiteral("Kr"), 83.798 },
            { QStringLiteral("Rb"), 85.4678 },
            { QStringLiteral("Sr"), 87.62 },
            { QStringLiteral("Y"),  88.90584 },
            { QStringLiteral("Zr"), 91.224 },
            { QStringLiteral("Nb"), 92.90637 },
            { QStringLiteral("Mo"), 95.95 },
            { QStringLiteral("Tc"), 98.0 },
            { QStringLiteral("Ru"), 101.07 },
            { QStringLiteral("Rh"), 102.90550 },
            { QStringLiteral("Pd"), 106.42 },
            { QStringLiteral("Ag"), 107.8682 },
            { QStringLiteral("Cd"), 112.414 },
            { QStringLiteral("In"), 114.818 },
            { QStringLiteral("Sn"), 118.710 },
            { QStringLiteral("Sb"), 121.760 },
            { QStringLiteral("Te"), 127.60 },
            { QStringLiteral("I"),  126.90447 },
            { QStringLiteral("Xe"), 131.293 },
            { QStringLiteral("Cs"), 132.90545 },
            { QStringLiteral("Ba"), 137.327 },
            { QStringLiteral("La"), 138.90547 },
            { QStringLiteral("Ce"), 140.116 },
            { QStringLiteral("Pr"), 140.90766 },
            { QStringLiteral("Nd"), 144.242 },
            { QStringLiteral("Pm"), 145.0 },
            { QStringLiteral("Sm"), 150.36 },
            { QStringLiteral("Eu"), 151.964 },
            { QStringLiteral("Gd"), 157.25 },
            { QStringLiteral("Tb"), 158.92535 },
            { QStringLiteral("Dy"), 162.500 },
            { QStringLiteral("Ho"), 164.93033 },
            { QStringLiteral("Er"), 167.259 },
            { QStringLiteral("Tm"), 168.93422 },
            { QStringLiteral("Yb"), 173.045 },
            { QStringLiteral("Lu"), 174.9668 },
            { QStringLiteral("Hf"), 178.49 },
            { QStringLiteral("Ta"), 180.94788 },
            { QStringLiteral("W"),  183.84 },
            { QStringLiteral("Re"), 186.207 },
            { QStringLiteral("Os"), 190.23 },
            { QStringLiteral("Ir"), 192.217 },
            { QStringLiteral("Pt"), 195.084 },
            { QStringLiteral("Au"), 196.966569 },
            { QStringLiteral("Hg"), 200.592 },
            { QStringLiteral("Tl"), 204.38 },
            { QStringLiteral("Pb"), 207.2 },
            { QStringLiteral("Bi"), 208.98040 },
            { QStringLiteral("Po"), 209.0 },
            { QStringLiteral("At"), 210.0 },
            { QStringLiteral("Rn"), 222.0 },
            { QStringLiteral("Fr"), 223.0 },
            { QStringLiteral("Ra"), 226.0 },
            { QStringLiteral("Ac"), 227.0 },
            { QStringLiteral("Th"), 232.0377 },
            { QStringLiteral("Pa"), 231.03588 },
            { QStringLiteral("U"),  238.02891 },
            { QStringLiteral("Np"), 237.0 },
            { QStringLiteral("Pu"), 244.0 },
            { QStringLiteral("Am"), 243.0 },
            { QStringLiteral("Cm"), 247.0 },
            { QStringLiteral("Bk"), 247.0 },
            { QStringLiteral("Cf"), 251.0 },
            { QStringLiteral("Es"), 252.0 },
            { QStringLiteral("Fm"), 257.0 },
            { QStringLiteral("Md"), 258.0 },
            { QStringLiteral("No"), 259.0 },
            { QStringLiteral("Lr"), 262.0 },
            { QStringLiteral("Rf"), 267.0 },
            { QStringLiteral("Db"), 268.0 },
            { QStringLiteral("Sg"), 269.0 },
            { QStringLiteral("Bh"), 270.0 },
            { QStringLiteral("Hs"), 269.0 },
            { QStringLiteral("Mt"), 278.0 },
            { QStringLiteral("Ds"), 281.0 },
            { QStringLiteral("Rg"), 282.0 },
            { QStringLiteral("Cn"), 285.0 },
            { QStringLiteral("Nh"), 286.0 },
            { QStringLiteral("Fl"), 289.0 },
            { QStringLiteral("Mc"), 290.0 },
            { QStringLiteral("Lv"), 293.0 },
            { QStringLiteral("Ts"), 294.0 },
            { QStringLiteral("Og"), 294.0 }
        };

        return masses;
    }
};

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QWidget window;
    window.setWindowTitle(QStringLiteral("Molar Mass Calculator"));
    window.setMinimumWidth(520);

    QLabel *titleLabel = new QLabel(QStringLiteral("Chemical formula:"));
    QLineEdit *formulaEdit = new QLineEdit;
    formulaEdit->setPlaceholderText(QStringLiteral("Examples: H2O, Ca(OH)2, CuSO4.5H2O"));
    formulaEdit->setFocus();

    QPushButton *calcButton = new QPushButton(QStringLiteral("Calculate"));
    calcButton->setDefault(true);

    QLabel *resultLabel = new QLabel(QStringLiteral("Enter a formula and press Calculate."));
    resultLabel->setWordWrap(true);

    QVBoxLayout *layout = new QVBoxLayout(&window);
    layout->addWidget(titleLabel);
    layout->addWidget(formulaEdit);
    layout->addWidget(calcButton);
    layout->addWidget(resultLabel);

    QObject::connect(calcButton, &QPushButton::clicked, [formulaEdit, resultLabel]() {
        QString error;
        double mass = 0.0;

        if (FormulaParser::parse(formulaEdit->text(), mass, error)) {
            QString massText = QString::number(mass, 'f', 6);

            // Remove useless trailing zeros, e.g. 18.015000 -> 18.015
            while (massText.endsWith(QChar('0')))
                massText.chop(1);

            if (massText.endsWith(QChar('.')))
                massText.chop(1);

            resultLabel->setText(QStringLiteral("%1 g/mol").arg(massText));
        } else {
            resultLabel->setText(QStringLiteral("Error: %1").arg(error));
        }
    });

    QObject::connect(formulaEdit,
                     &QLineEdit::returnPressed,
                     calcButton,
                     &QPushButton::click);

    window.show();

    return app.exec();
}
