#include <BedRegion.h>

BedRegion::BedRegion() {
  chr = "chr0";
  start = 0;
  end = 0;
}

BedRegion::BedRegion(std::string _chr, int _start, int _end) {
  chr = _chr;
  start = _start;
  end = _end;
}

bool BedRegion::parse(std::string str) {
  char chr[40];
  int st = -1, end = -1;
  // Try chr:start-end format first
  sscanf(str.c_str(), "%30[^:]:%d-%d", chr, &st, &end);
  // If that fails, try chr:start:end format (colon separator)
  if (st == -1 || end == -1) {
    st = -1; end = -1;
    sscanf(str.c_str(), "%30[^:]:%d:%d", chr, &st, &end);
  }
  if (st != -1 && end != -1) {
    this->start = st;
    this->end = end;
    this->chr = (std::string)chr;
    return true;
  }

  return false;
}

void BedRegion::print() { printf("%s %d %d\n", chr.c_str(), start, end); }

bool BedRegion::contains(int pos) { return (pos >= start && pos <= end); }

bool BedRegion::tryParse(std::string str) {
  char chr[40];
  int st = -1, en = -1;
  // Try chr:start-end format first
  sscanf(str.c_str(), "%30[^:]:%d-%d", chr, &st, &en);
  // If that fails, try chr:start:end format (colon separator)
  if (st == -1 || en == -1) {
    st = -1; en = -1;
    sscanf(str.c_str(), "%30[^:]:%d:%d", chr, &st, &en);
  }
  if (st == -1 || en == -1)
    return false;
  return true;
}
