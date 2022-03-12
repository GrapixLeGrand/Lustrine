#include "JsonWriter.h"

JsonWriter &JsonWriter::item() {
    if (first) first = false;
    else os << ',';
    return *this;
}

JsonWriter::JsonWriter(std::ostream &os) : os{os} {}

JsonWriter &JsonWriter::value(float v) {
    item();
    os << v;
    return *this;
}
JsonWriter &JsonWriter::value(bool v) {
    item();
    os << (v ? "true" : "false");
    return *this;
}
JsonWriter &JsonWriter::value(uint32_t v) {
    item();
    os << v;
    return *this;
}
JsonWriter &JsonWriter::value(int32_t v) {
    item();
    os << v;
    return *this;
}
JsonWriter &JsonWriter::value(const char *v) {
    item();
    os << '"';
    for (const char *c = v; *c != '\0'; ++c) {
        if (*c == '\\') os << "\\\\";
        else if (*c == '"') os << "\\\"";
        else if (*c < ' ') {
            if (*c == '\n') os << "\\n";
            else if (*c == '\r') os << "\\r";
            else if (*c == '\t') os << "\\t";
        } else if (*c == '\x7F');
        else os << *c;
    }
    os << '"';
    return *this;
}
JsonWriter &JsonWriter::null() {
    item();
    os << "null";
    return *this;
}
JsonWriter &JsonWriter::safeKey(const char *v) {
    item();
    os << '"' << v << "\":";
    return *this;
}
JsonWriter &JsonWriter::unsafeKey(const char *v) {
    item();
    value(v);
    os << ':';
    return *this;
}
JsonWriter &JsonWriter::bObject() {
    item();
    os << '{';
    first = true;
    return *this;
}
JsonWriter &JsonWriter::eObject() {
    os << '}';
    return *this;
}
JsonWriter &JsonWriter::bArray() {
    item();
    os << '[';
    first = true;
    return *this;
}
JsonWriter &JsonWriter::eArray() {
    os << ']';
    return *this;
}