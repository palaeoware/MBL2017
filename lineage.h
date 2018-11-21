#ifndef LINEAGE_H
#define LINEAGE_H

#include <QObject>
#include <QString>
#include <QLinkedList>
#include "genus.h"

//#define CHARACTER_WORDS 4 //number of 32 bit words used for characters. Default 4 gives 128 characters
//#define CHARACTER_WORDS 8 //256 characters
//#define CHARACTER_WORDS 16 //512 characters
#define CHARACTER_WORDS 64 //1024 characters
//memory usage during a run is proportional to CHARACTER_WORDS

class Lineage
{
public:
    qint64 id;
    int simple_id; //used for some tree dumps that want numbers 1-n with no gaps for extant species
    Lineage(quint32 characters[], Lineage *parent, qint64 timestamp, quint32 initial_characters[]=(quint32*)0);
    ~Lineage();
    qint64 time_created;
    qint64 time_died;
    qint64 time_split;
    qint64 genusnumber;
    Lineage *daughter_lineage_A;
    Lineage *daughter_lineage_B;
    Lineage *parent_lineage;
    bool dontpropogatedelete;
    quint32 initial_characters[CHARACTER_WORDS];
    quint32 current_characters[CHARACTER_WORDS];


    QString dump();
    void iterate(qint64 timestamp);
    int count_alive();
    int count_extinct();
    int count_branched();

    QString newickstring();
    void genuslabels(qint64 currentlabel);
    QString xmlstring();
    void getextantlist(QList<Lineage *> *list);
    void genuslabels_distance(qint64 currentlabel, int treelength);
    QString tntstring();
    bool isThisADescendent(Lineage *desc);
    Lineage *strip_extinct(Lineage *parent);
    void cull_dead_branches();
    Lineage *getsister();
    bool RDT_check(quint64 currenttime);
    void RDT_incorporate(Genus *g);
    void resetGenusLabels();
    void genuslabels_IDT(qint64 currentgenusnumber, int currenttime);
    void getleaflist(QLinkedList<Lineage *> *list);
    void addtogenus(qint64 genusnum);
    int nodelength(int currenttime);
    bool DoesThisLeafSetContainACladeDescendingFromMe(QSet<Lineage *> *potentialspecieslist);
    bool isThisAnAncestor(Lineage *desc);
    QString numbertree(int add);
    void dosimplenumbers();
    Lineage *find_clade_with_precise_size(int preciseleafcount);
    QString getcharactermatrix(bool root);
    void addstringgenomestoset(QSet<QString> *set);
    int count_ones(int word, quint32 mask);
    int count_zeros(int word, quint32 mask);
    QString getcharactersasstring(quint32 characters[]);
private:

    void setgenuslabels(qint64 currentlabel, bool labelfossils=false);
    int persistedfor();
    void domutation();
};

extern int CHANCE_EXTINCT;
extern int CHANCE_SPECIATE;
extern double CHANCE_EXTINCT_DOUBLE;
extern double CHANCE_SPECIATE_DOUBLE;
extern quint32 CHANCE_MUTATE;
extern int ABS_THRESHOLD;
extern qint64 speccount;
extern qint64 diedcount;
extern qint64 totaltries;
#endif // LINEAGE_H
