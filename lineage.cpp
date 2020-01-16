#include "lineage.h"
#include <QTextStream>
#include "simulation.h"
#include <QDebug>
#include <QList>
#include "genus.h"


/////////////////////////////////////////////////////////////////////
//Lineage class - equates essentially to a species. Lineages have
//a timestamp for when they appeared, went extinct, or split
//also holds pointers to parent and two daughter lineages, initialised to 0
//timestamps are initialised to -1, and remain so if unused (e.g. if lineage
//goes extinct and does not split, time_split remains -1. Split lineages
//are replaced by their daughters in the simulation
//dontpropogatedelete is a flag used in stripping extinct genera to stop daughter lineages being
//deleted before they can be reassigned.
/////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////
//Imports of globals from simulations.cpp
/////////////////////////////////////////////////////
extern double RDT_THRESHOLD;
extern int IDT_THRESHOLD;
extern qint64 genusnumberidt;
extern Lineage *dummy_parameter_lineage;


/////////////////////////////////////////////////////
//Constructor/destructor
/////////////////////////////////////////////////////

Lineage::Lineage(quint32 characters[], Lineage *parent, qint64 timestamp, quint32 initialchars[])
//initial characters is default 0 pointer, in which case none supplied, use current
{

    //tracked=false;
    //if (parent==0) tracked=true;
    time_created=timestamp;
    time_died=-1;
    time_split=-1;
    daughter_lineage_A=nullptr;
    daughter_lineage_B=nullptr;
    genusnumber=0;
    parent_lineage=parent;
    simple_id=-1;

    dontpropogatedelete=false;
    id=TheSimGlobal->nextid++;

    //copy characters
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        current_characters[i]=characters[i];
        if (initialchars)
            initial_characters[i]=initialchars[i];
        else
            initial_characters[i]=characters[i];
    }

}

Lineage::~Lineage()
{
    if (dontpropogatedelete) return; //switch to turn off deletion of children
    //Recursively delete all daughter lineages via their destructors;
    if (daughter_lineage_A!=nullptr) delete daughter_lineage_A;
    if (daughter_lineage_B!=nullptr) delete daughter_lineage_B;
}

/////////////////////////////////////////////////////
//Iteration - MBL algorithm implementation
/////////////////////////////////////////////////////

void Lineage::domutation()
{
    for (int i=0; i<CHARACTER_WORDS*32; i++)
    {
        if (TheSimGlobal->GetRandom()<CHANCE_MUTATE)
            current_characters[i/32] ^= tweakers[i%32];
    }
}

void Lineage::iterate(qint64 timestamp)
{
    //iterate a lineage - i.e. perform MBL calculations
    if (time_died!=-1) return; //it's dead - no need to do anything
    if (time_split!=-1) //it has split - iterate it's daughters instead
    {
        daughter_lineage_A->iterate(timestamp);
        daughter_lineage_B->iterate(timestamp);
        return;
    }

    domutation();
    //get a random double (0-1)
    quint32 r=TheSimGlobal->GetRandom();
    double rdouble=((double)r)/65536;
    rdouble/=65536;

    if (rdouble<CHANCE_EXTINCT_DOUBLE) //lineage went extinct
    {
        TheSimGlobal->leafcount--;
        time_died=timestamp; //record when this happened
        return;
    }

    if (rdouble<(CHANCE_EXTINCT_DOUBLE+CHANCE_SPECIATE_DOUBLE)) //lineage speciated
    {
        TheSimGlobal->leafcount++;
        time_split=timestamp; //record when this happened
        daughter_lineage_A=new Lineage(current_characters, this,timestamp); //create two new daughter lineages with same characters
        daughter_lineage_B=new Lineage(current_characters, this,timestamp);
        return;
    }
}


/////////////////////////////////////////////////////
//Recursive functions to report on/summarise tree
/////////////////////////////////////////////////////

int Lineage::count_alive()
{
    //recursively count how many descendents of this lineage are alive (don't have a time_died set) - including this one
    if (daughter_lineage_A && daughter_lineage_B) //it's split - value is count_alive of both branches
        return daughter_lineage_A->count_alive()+daughter_lineage_B->count_alive();

    if (time_died!=-1) return 0;
    return 1;
}

int Lineage::count_extinct()
{
    //recursively count how many descendents of this lineage are extinct (DO have a time_died set) - including this one
    if (daughter_lineage_A && daughter_lineage_B) //it's split
        return daughter_lineage_A->count_extinct()+daughter_lineage_B->count_extinct();

    if (time_died==-1) return 0; //should be impossible anyway - it can't be split and extinct
    return 1;
}

int Lineage::count_branched()
{
    //recursively count how many descendents of this lineage speciated - including this one
    if (daughter_lineage_A && daughter_lineage_B) //it's split
        return 1+daughter_lineage_A->count_branched()+daughter_lineage_B->count_branched();

    return 0;
}

void Lineage::getextantlist(QList<Lineage *> *list)
{
    //appends all descending extant species to the list (including this one_
    if (time_split==-1 && time_died==-1) //extant
    {
        list->append(this);
    }
    else
    {
        if (time_split!=-1) //if it split - recurse onto daughters
        {
            daughter_lineage_A->getextantlist(list);
            daughter_lineage_B->getextantlist(list);
        }
    }
}

void Lineage::getleaflist(QLinkedList<Lineage *> *list) //includes extinct
{
    //appends all descending extant or extinct species to the list (including this one)
    if (time_split==-1) //a leaf
    {
        list->append(this);
    }
    else
    {
        daughter_lineage_A->getleaflist(list);
        daughter_lineage_B->getleaflist(list);
    }
}

bool Lineage::isThisADescendent(Lineage *desc)
{
    //recurseive function - does the lineage passed appear in descendents of this lineage?
    if (this==desc) return true;
    if (time_split==-1) return false; //no descendents and it wasn't this, so false
    else
    {
        bool b1=daughter_lineage_A->isThisADescendent(desc);
        bool b2=daughter_lineage_B->isThisADescendent(desc);
        if (b2 || b1) return true; else return false;
    }
}

bool Lineage::isThisAnAncestor(Lineage *desc)
{
    //recursive function - does the lineage passed appear in ancestors of this lineage?
    if (this==desc) return true;
    if (!parent_lineage) return false; //hit root - didn't find it
    return parent_lineage->isThisAnAncestor(desc);
}


bool Lineage::DoesThisLeafSetContainACladeDescendingFromMe(QSet<Lineage *> *potentialspecieslist)
{
    //er... you can guess this one
    if (time_split==-1)
    {
        if (potentialspecieslist->contains(this)) return true; //fine, I am in the clade
        else return false; //ah, no I'm not in the list, so I am an excluded species
    }
    else
    {
        bool b1=daughter_lineage_A->DoesThisLeafSetContainACladeDescendingFromMe(potentialspecieslist);
        bool b2=daughter_lineage_B->DoesThisLeafSetContainACladeDescendingFromMe(potentialspecieslist);
        if (b2 && b1) return true; else return false;
    }
}

void Lineage::resetGenusLabels()
{
    genusnumber=0;
    if (time_split!=-1)
    {
        daughter_lineage_A->resetGenusLabels();
        daughter_lineage_B->resetGenusLabels();
    }
}

/////////////////////////////////////////////////////
//Generate tree-file output recursively
/////////////////////////////////////////////////////


QString Lineage::tntstring()
{
    //recursively generate TNT-format text description of tree
    if (time_split==-1)
    {
        QString s;
        if (time_died==-1)
            s.sprintf("%lld",id);
        return s;
    }
    else
    {
        QString s;
        QTextStream out(&s);
        int A=daughter_lineage_A->count_alive();
        int B=daughter_lineage_B->count_alive();

        if (A==0 && B==0) s="";
        if (A>0 && B>0)
        {
            out<<"("<<daughter_lineage_A->tntstring()<<" "<<daughter_lineage_B->tntstring()<<")";
        }
        if (A>0 && B==0)
        {
            return daughter_lineage_A->tntstring();
        }
        if (B>0 && A==0)
        {
            return daughter_lineage_B->tntstring();
        }
        return s;
    }
}

QString Lineage::numbertree(int add)
{
    //recursively generate number-format text description of tree
    int bl=TheSimGlobal->currenttime-time_created;
    if (time_died!=-1) bl=time_died-time_created;
    if (time_split!=-1) bl=time_split-time_created;

    if (time_split==-1)
    {
        return QString("%1:%2").arg(simple_id+add).arg(bl);

    }
    else
    {
        QString s;
        QTextStream out(&s);
        out<<"("<<daughter_lineage_A->numbertree(add)<<","<<daughter_lineage_B->numbertree(add)<<"):"<<bl;
        return s;
    }
}

QString Lineage::newickstring()
{
    //recursively generate Newick-format text description of tree
    int bl=TheSimGlobal->currenttime-time_created;
    if (time_died!=-1) bl=time_died-time_created;
    if (time_split!=-1) bl=time_split-time_created;

    if (time_split==-1)
    {
        QString s;
        if (time_died==-1)
            s.sprintf("g%lld_s%lld:%d",genusnumber,id,bl);
        else
            s.sprintf("{g%lld_s%lld}:%d",genusnumber,id,bl);
        return s;
    }
    else
    {
        QString s;
        QTextStream out(&s);
        out<<"("<<daughter_lineage_A->newickstring()<<","<<daughter_lineage_B->newickstring()<<"):"<<bl;
        return s;
    }
}

QString Lineage::xmlstring()
{
    //recursively generate xml-format text description of tree
    int bl=TheSimGlobal->currenttime-time_created;
    if (time_died!=-1) bl=time_died-time_created;
    if (time_split!=-1) bl=time_split-time_created;

    if (time_split==-1)
    {
        QString s;
        if (time_died==-1)
            s.sprintf("<clade><name>g%lld_s%lld</name><branch_length>%d</branch_length></clade>\n",genusnumber,id,bl);
        else
            s.sprintf("<clade><name>Ext_s%lld</name><branch_length>%d</branch_length></clade>\n",id,bl);
        return s;
    }
    else
    {
        QString s;
        QTextStream out(&s);
        out<<"<clade><name>n"<<id<<"</name><branch_length>"<<bl<<"</branch_length>"<<daughter_lineage_A->xmlstring()<<daughter_lineage_B->xmlstring()<<"</clade>\n";
        return s;
    }
}


void Lineage::dosimplenumbers()
{
    //for tree algorithms really. Provides secondary numbering system from 1-n, with no gaps
    if (time_split==-1) //a leaf
        simple_id=TheSimGlobal->nextsimpleIDnumber++;
    else
    {
        daughter_lineage_A->dosimplenumbers();
        daughter_lineage_B->dosimplenumbers();
    }
}

/////////////////////////////////////////////////////
//IDT algorithm
/////////////////////////////////////////////////////

int Lineage::nodelength(int currenttime)
{
    //branchlength to present or split or extinction
    if (time_split==-1)
    {
        if (time_died==-1)
            return currenttime-time_created;
        else
            return time_died-time_created;
    }
    else
        return time_split-time_created;
}

void Lineage::genuslabels_IDT(qint64 currentgenusnumber, int currenttime)
{
    // - look at the two lengths for daughters. if EITHER is over threshold, both get diffent genus numbers
    // - if this is not split, just add to node and return as before

    if (time_split==-1) // a leaf - extinct or otherwise - just add to current genus
    {
        genusnumber=currentgenusnumber;
        Genus *g;
        if (!(genera.contains(currentgenusnumber))) //genus doesn't exist - create and add to genera
        {
            g=new Genus;
            g->id=currentgenusnumber;
            genera.insert(g->id,g);
        }
        else
            g=genera.value(currentgenusnumber); //already exists

        g->species.append(this); //add this species to genus
    }
    else // a split.
    {
        int d1=daughter_lineage_A->nodelength(currenttime);
        int d2=daughter_lineage_B->nodelength(currenttime);

        if (d1>IDT_THRESHOLD || d2>IDT_THRESHOLD)
        {
            currentgenusnumber=++genusnumberidt;
            daughter_lineage_A->genuslabels_IDT(currentgenusnumber,currenttime);
            currentgenusnumber=++genusnumberidt;
            daughter_lineage_B->genuslabels_IDT(currentgenusnumber,currenttime);
        }
        else
        {
            //node was not too long - keep same genus number for both children
            daughter_lineage_A->genuslabels_IDT(currentgenusnumber,currenttime);
            daughter_lineage_B->genuslabels_IDT(currentgenusnumber,currenttime);
        }
    }
}




/////////////////////////////////////////////////////
//FDT algorithm
/////////////////////////////////////////////////////

void Lineage::genuslabels_distance(qint64 currentlabel, int treelength)
{
    //used to apply FDT to whole tree. Recursive, normally called on rootspecies first
    if (time_split==-1)
    {
        //Not split, so extant or extinct

        genusnumber=currentlabel; //put in a genus - applies also to extinct, though we do nothing with those

        if (time_died==-1) //extant
        {
            Genus *g;
            if (!(genera.contains(currentlabel))) //genus doesn't exist - create and add to genera
            {
                g=new Genus;
                g->id=currentlabel;
                genera.insert(g->id,g);
            }
            else
                g=genera.value(currentlabel); //already exists

            g->species.append(this); //add this species to genus
        }
        //do nothing if extinct
        return;
    }
    else
    {
        //it split - see if split is before or after the threshold
        if ((treelength-time_split) > ABS_THRESHOLD)
        {
            //before - so the daughters are in separate genera
            daughter_lineage_A->genuslabels_distance(TheSimGlobal->nextgenusnumber++,treelength);
            daughter_lineage_B->genuslabels_distance(TheSimGlobal->nextgenusnumber++,treelength);
        }
        else
        {
            //after - so in same genera. Set labels of all children
            daughter_lineage_A->setgenuslabels(currentlabel);
            daughter_lineage_B->setgenuslabels(currentlabel);
        }
    }
}

void Lineage::setgenuslabels(qint64 currentlabel, bool labelfossils)
{
    //recursively set all species descending from here to the passed genus label
    genusnumber=currentlabel;
    if (time_split!=-1)
    {
        daughter_lineage_A->setgenuslabels(currentlabel,labelfossils);
        daughter_lineage_B->setgenuslabels(currentlabel,labelfossils);
    }
    else
    {
        if (time_died==-1 || labelfossils) //if labelfossils on - well label it even if it did die!
        {
            Genus *g;
            if (!(genera.contains(currentlabel)))
            {
                g=new Genus;
                g->id=currentlabel;
                genera.insert(g->id,g);
            }
            else
                g=genera.value(currentlabel);

            g->species.append(this);
        }
    }
}


/////////////////////////////////////////////////////
//Recursive functions to remove extinct taxa from tree prior to RDT. Returns pointer to new root.
/////////////////////////////////////////////////////

void Lineage::cull_dead_branches()
{
    //recurse through tree. Where we find a branch time but only a single branch - merge the two structures, deleting the upper one
    if (time_split!=-1) //do nothing if not a split
    {
        if (daughter_lineage_A==nullptr && daughter_lineage_B==nullptr)
        {
            qDebug()<<"Internal error in cull_dead_branches. Oh dear.";
        }
        if (daughter_lineage_A==nullptr) //must mean B was not 0 - both 0 case was culled in strip_extinct
        {
            time_died=daughter_lineage_B->time_died;
            time_split=daughter_lineage_B->time_split;
            Lineage *todelete=daughter_lineage_B;
            daughter_lineage_A=todelete->daughter_lineage_A;
            daughter_lineage_B=todelete->daughter_lineage_B;
            if (daughter_lineage_A) daughter_lineage_A->parent_lineage=this;
            if (daughter_lineage_B) daughter_lineage_B->parent_lineage=this;

            todelete->dontpropogatedelete=true;
            delete todelete;
            //parent and time created stay the same
            cull_dead_branches(); //do this again on me
            return;
        }
        if (daughter_lineage_B==nullptr) //must mean A was not 0 - both 0 case was culled in strip_extinct
        {
            time_died=daughter_lineage_A->time_died;
            time_split=daughter_lineage_A->time_split;
            Lineage *todelete=daughter_lineage_A;
            daughter_lineage_A=todelete->daughter_lineage_A;
            daughter_lineage_B=todelete->daughter_lineage_B;
            todelete->dontpropogatedelete=true;
            if (daughter_lineage_A) daughter_lineage_A->parent_lineage=this;
            if (daughter_lineage_B) daughter_lineage_B->parent_lineage=this;
            delete todelete;
            //parent and time created stay the same
            cull_dead_branches(); //do this again on me
            return;
        }
        //if here - it's a normal split
        daughter_lineage_A->cull_dead_branches();
        daughter_lineage_B->cull_dead_branches();
    }
}

Lineage *Lineage::strip_extinct(Lineage *parent)
{
    //return 0 pointer if there are no extant descendants
    //otherwise return pointer to new Lineage structure with extant only descs
    if (time_died==-1 && time_split==-1) //extant, didn't split
    {
        //create new lineage for this, return;
        Lineage *l = new Lineage(current_characters,parent,time_created,initial_characters);
        l->id=id;
        //everything else can stay at default
        return l;
    }
    if (time_died!=-1) //extinct branch
    {
        return (Lineage *)nullptr; //return null pointer - this branch is not real
    }
    if (time_split!=-1) //split node
    {
        Lineage *l = new Lineage(current_characters,parent,time_created,initial_characters);
        l->id=id;

        l->time_died=-1;
        l->time_split=time_split;
        l->daughter_lineage_A=daughter_lineage_A->strip_extinct(l);
        l->daughter_lineage_B=daughter_lineage_B->strip_extinct(l);
        //if both returned 0 - no extant descendants - we return 0 too
        if (l->daughter_lineage_A==nullptr && l->daughter_lineage_B==nullptr)
        {
            delete l;
            return (Lineage *)nullptr;
        }
//        if (l->daughter_lineage_A) l->daughter_lineage_A->parent_lineage=this;
//        if (l->daughter_lineage_B) l->daughter_lineage_B->parent_lineage=this;

        return l;
    }
    qDebug()<<"Error in strip extinct - unhandled case";
    return (Lineage *)nullptr;

}


/////////////////////////////////////////////////////
//Functions used by RDT algoritm
/////////////////////////////////////////////////////

Lineage * Lineage::getsister()
{
    //determine the sister clade to this clade - returns 0 if no sister clade (root)
    if (parent_lineage==nullptr) return nullptr;

    //must be one of parents daughters - find which
    if (parent_lineage->daughter_lineage_A==this)
        return parent_lineage->daughter_lineage_B;
    else
        return parent_lineage->daughter_lineage_A;
}

bool Lineage::RDT_check(quint64 currenttime)
{
    //checks to see whether this lineage could be incorporated into a genus
    //This function is documented in the manuscript
    if (time_split==-1) //already culled fossils, so must be a leaf
        return true;
    else
    {
        if ((quint64)(((double)(currenttime-time_split))/RDT_THRESHOLD)>=currenttime-time_created)
        {
            if (daughter_lineage_A->RDT_check(currenttime) && daughter_lineage_B->RDT_check(currenttime))
                return true;
            else
                return false;
        }
        else return false; //fails 50% age rule
    }
}

void Lineage::RDT_incorporate(Genus *g)
{
    //incorporate all tips descending from here into this genus
    if (time_split==-1)
    {
        g->species.append(this);
        genusnumber=g->id;
    }
    else
    {
        daughter_lineage_A->RDT_incorporate(g);
        daughter_lineage_B->RDT_incorporate(g);
    }
}

/////////////////////////////////////////////////////
//Find clade of a certain size within the tree
/////////////////////////////////////////////////////

Lineage * Lineage::find_clade_with_precise_size(int preciseleafcount)
{
    if (time_split==-1) return nullptr; //no daughters
    int a=daughter_lineage_A->count_alive();
    if (a==preciseleafcount) return daughter_lineage_A;
    int b=daughter_lineage_A->count_alive();
    if (b==preciseleafcount) return daughter_lineage_B;
    if (a>preciseleafcount)
    {
        Lineage *r=daughter_lineage_A->find_clade_with_precise_size(preciseleafcount);
        if (r!=nullptr) return r;
    }
    if (b>preciseleafcount)
    {
        Lineage *r=daughter_lineage_B->find_clade_with_precise_size(preciseleafcount);
        if (r!=nullptr) return r;
    }

    //if it gets here it's not found big enough daughter clades
    return nullptr;
}



/////////////////////////////////////////////////////
//Debugging or reporting functions
/////////////////////////////////////////////////////

int Lineage::persistedfor()
{
    //returns number of steps the lineage persisted for until
    //simulation ended or it split or it died
    if (time_died!=-1)
        return time_died-time_created;
    if (time_split!=-1)
        return time_split-time_created;

    //still alive - can't report this, don't know total number of iterations here!
    return -1;

}

QString Lineage::dump()
{
    QString output;
    QTextStream out(&output);
    out<<"<br /><br />Lineage dump<br />";
    out<<"timstamp created: "<<time_created<<"<br />";
    if (parent_lineage==nullptr) out<<"No parent lineage - root species<br />";
    else
        out<<"Parent lineage pointer: "<<parent_lineage<<"\n";

    if (time_died==-1 && time_split==-1) out<<"Status: Alive"<<"<br />";
    if (time_died!=-1) out<<"Status: Extinct, died at "<<time_died<<"<br />";
    if (time_split!=-1) out<<"Status: Split at "<<time_split<<")"<<"<br />";
    out<<"initial characters<br />";
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        QString binary;
        binary=QString::number(initial_characters[i],2);
        while (binary.length()<32)
            binary.prepend("0");
        out<<binary;
    }
    out<<"<br />";
    out<<"current characters<br />";
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        QString binary;
        binary=QString::number(current_characters[i],2);
        while (binary.length()<32)
            binary.prepend("0");
        out<<binary;
    }
    out<<"<br />";

    return output;
}
QString Lineage::getcharactersasstring(quint32 characters[])
{
    QString output;
    QTextStream out(&output);
    for (int i=0; i<CHARACTER_WORDS; i++)
    {
        QString binary;
        binary=QString::number(characters[i],2);
        while (binary.length()<32)
            binary.prepend("0");
        out<<binary;
    }
    return output;
}

int Lineage::count_zeros(int word, quint32 mask)
{

    if (time_split!=-1) //split - recurse
    {
        return daughter_lineage_A->count_zeros(word,mask)+daughter_lineage_B->count_zeros(word,mask);
    }
    else
    {
        if ((current_characters[word]&mask)==0) return 1; else return 0;
    }
}

int Lineage::count_ones(int word, quint32 mask)
{

    if (time_split!=-1) //split - recurse
    {
        return daughter_lineage_A->count_ones(word,mask)+daughter_lineage_B->count_ones(word,mask);
    }
    else
    {
        if ((current_characters[word]&mask)==1) return 1; else return 0;
    }
}


void Lineage::addstringgenomestoset(QSet<QString> *set)
{
    if (time_split!=-1) //split - recurse
    {
        daughter_lineage_A->addstringgenomestoset(set);
        daughter_lineage_B->addstringgenomestoset(set);
    }
    else
    set->insert(getcharactersasstring(current_characters));
}


QString Lineage::getcharactermatrix(bool root)
{
    QString retval;

    if (root)
    //this is root species - so outgroup - need starting characters
    {
        retval= QString("Species_%1        ").arg(0).left(14);
        retval+=getcharactersasstring(initial_characters);
        retval+="\n";
    }
    if (time_split!=-1) //split - recurse
    {
        return retval+daughter_lineage_A->getcharactermatrix(false)+daughter_lineage_B->getcharactermatrix(false);
    }
    return retval+QString("Species_%1        ").arg(simple_id).left(14)+getcharactersasstring(current_characters)+"\n";
}

