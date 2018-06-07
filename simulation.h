#ifndef SIMULATION_H
#define SIMULATION_H
#include <QObject>
#include <QHash>
#include <QVector>

#define TREE_MODE_UNCLASSIFIED 0
#define TREE_MODE_FDT 1
#define TREE_MODE_RDT 2
#define TREE_MODE_IDT 3
#define TREE_MODE_TNTMB 9 //just for the TNT/MrBayes mode - internal only


//Define below for TREE_MODE_MAX should be the max value used above
#define TREE_MODE_MAX 8

#define PROPORTIONAL_BINS 20 //100/this must be integral. So 10, 20, 25, 50 probably. Actually have this+1 bins - there is one for 100%

class Lineage;
class MainWindow;
class Genus;
class Simulation
{
public:
    Simulation();
    ~Simulation();
    void run(MainWindow *mainwin);
    MainWindow *mw;
    Lineage *rootspecies;
    Lineage *crownroot; //root for extant-only version of tree
    qint64 currenttime;
    void increment();
    quint32 randoms[65536];
    quint16 rpoint;
    quint64 rfilepoint;
    quint32 GetRandom();
    qint64 nextid;
    qint64 nextgenusnumber;
    bool stopflag;
    int leafcount;
    int nextsimpleIDnumber;

    void read65536numbers();
    quint32 GetRandom16();
    void stop();
    QString filepath;
    QString dumpnewick(Lineage *rootlineage);
    QString dumpphyloxml(Lineage *rootlineage);
    QString dumptnt(Lineage *rootlineage);
    QString modeToString(int mode);
    QList<QHash<int,int>*> counts;
    QList<QList<int>*> proportional_counts; //in PROPORTIONAL_BINS bins - 100/PROPORTIONAL_BINS must be an integer
    QList<QList<int>*> saturated_tree_sizes;
    QString modeToShortString(int mode);
    int distancebetween(quint32 chars1[], quint32 chars2[]);
    QString dumpnex_with_matrix(Lineage *rootlineage, int iter);
    QString dumptnt_with_matrix(Lineage *rootlineage, int iter);
    QString dump_nex_alone(Lineage *rootlineage);
    void randomcharacters(quint32 *chars);
private:
    quint32 rand256();

    int genera_data_report(int mode);
    void RDT_genera();
    Lineage *get_mrca(Lineage *l0, Lineage *l1);
    int get_mrca_age(QList<Lineage *> *terminals);
};

struct maxgenusdatapoint
{
    quint32 maxgenussize;
    quint32 treesize;
};

extern Simulation *TheSimGlobal;
extern quint32 tweakers[32];
extern int bitcounts[65536];
extern QHash<qint64,Genus*> genera;

#endif // SIMULATION_H
