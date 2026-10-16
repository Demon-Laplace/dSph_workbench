/* One-pass equivalent of two ICmerger_matrixMW calls.
 * Compatibility is intentional: float operations, 0.017453 rad/degree,
 * Rx/Ry/Rz order, ignored first-input angles, %.7g intermediate rounding,
 * stable type ordering and ID offsets. Do NOT enable fast-math/FMA.
 * The reference executable SHA256 and verification are documented in MERGE_IC.md.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

typedef struct { int id; float v[8]; int type; } Particle;
typedef struct { char input[2][4096]; float v[12]; } Parameters;
typedef struct { Particle *p; int counts[4]; int n; } Input;

static void fail(const char *message) { fprintf(stderr, "%s\n", message); exit(2); }
static FILE *open_read(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(2); }
    return f;
}
static Parameters parameters(const char *path) {
    Parameters a;
    FILE *f = open_read(path);
    if (fscanf(f, "%4095s %4095s", a.input[0], a.input[1]) != 2) fail("Invalid merger inputs");
    for (int j=0; j<12; j++)
        if (fscanf(f, "%f", &a.v[j]) != 1 || !isfinite(a.v[j])) fail("Invalid merger transform");
    fclose(f);
    return a;
}
static void matrices(const Parameters *a, float m[3][3][3]) {
    memset(m, 0, 27*sizeof(float));
    for (int axis=0; axis<3; axis++) {
        double angle = (double)a->v[9+axis] * 0.017453;
        float s = (float)sin(angle), c = (float)cos(angle);
        int i=(axis+1)%3, j=(axis+2)%3;
        m[axis][axis][axis]=1;
        m[axis][i][i]=m[axis][j][j]=c;
        m[axis][i][j]=s; m[axis][j][i]=-s;
    }
}
static void transform(Particle *p, const Parameters *a, float m[3][3][3]) {
    for (int offset=0; offset<=3; offset+=3) {
        for (int axis=0; axis<3; axis++) {
            float out[3]={0,0,0};
            for (int i=0; i<3; i++)
                for (int j=0; j<3; j++) out[i] += p->v[offset+j]*m[axis][i][j];
            memcpy(p->v+offset, out, sizeof(out));
        }
    }
    for (int j=0; j<6; j++) p->v[j] += a->v[j];
}
static void round_intermediate(Particle *p) {
    char buffer[64];
    for (int j=0; j<8; j++) {
        snprintf(buffer, sizeof(buffer), "%.7g", (double)p->v[j]);
        p->v[j]=strtof(buffer, NULL);
    }
}
static Input read_input(const char *path, int role, const Parameters *galaxy,
                        const Parameters *cube, float gm[3][3][3], float cm[3][3][3]) {
    Input a = {0};
    FILE *f=open_read(path);
    int64_t total=0;
    for (int j=0; j<4; j++) {
        if (fscanf(f, "%d", &a.counts[j]) != 1 || a.counts[j]<0) fail("Invalid particle counts");
        total+=a.counts[j];
    }
    if (total>INT_MAX) fail("Particle count exceeds legacy int32 range");
    a.n=(int)total;
    a.p=malloc((size_t)(a.n ? a.n : 1)*sizeof(Particle));
    if (!a.p) fail("Insufficient memory for particles");
    int seen[4]={0};
    for (int i=0; i<a.n; i++) {
        Particle *p=&a.p[i];
        if (fscanf(f, "%d %f %f %f %f %f %f %f %f %d", &p->id,
                   &p->v[0], &p->v[1], &p->v[2], &p->v[3], &p->v[4],
                   &p->v[5], &p->v[6], &p->v[7], &p->type) != 10) fail("Truncated/invalid particle row");
        if (p->id<0 || p->type<0 || p->type>3) fail("Unsupported particle ID/type");
        seen[p->type]++;
        for (int j=0; j<8; j++) if (!isfinite(p->v[j])) fail("Non-finite particle field");
        if (p->v[6]<=0) fail("Non-positive particle mass");
        if (role==2) transform(p, galaxy, gm);
        if (role!=0) { round_intermediate(p); transform(p, cube, cm); }
        for (int j=0; j<8; j++) if (!isfinite(p->v[j])) fail("Non-finite transformed field");
    }
    for (int j=0; j<4; j++) if (seen[j]!=a.counts[j]) fail("Particle type counts disagree with header");
    int ch;
    while ((ch=fgetc(f))!=EOF) if (!isspace(ch)) fail("Unexpected trailing particle data");
    fclose(f);
    return a;
}
static int compare_ids(const void *a, const void *b) {
    int x=*(const int *)a, y=*(const int *)b;
    return (x>y)-(x<y);
}
int main(int argc, char **argv) {
    if (argc!=4) fail("Usage: ic_merge mergegalaxy.par mergecube.par output.dat (paths relative to cwd)");
    Parameters g=parameters(argv[1]), c=parameters(argv[2]);
    float gm[3][3][3], cm[3][3][3];
    matrices(&g,gm); matrices(&c,cm);
    /* Final ordering: type, then cube / first galaxy / second galaxy. */
    Input a[3];
    a[0]=read_input(c.input[0],0,&g,&c,gm,cm);
    a[1]=read_input(g.input[0],1,&g,&c,gm,cm);
    a[2]=read_input(g.input[1],2,&g,&c,gm,cm);
    int64_t total=(int64_t)a[0].n+a[1].n+a[2].n;
    if (total<=0 || total>INT_MAX) fail("Invalid combined particle count");
    int *ids=malloc((size_t)total*sizeof(int));
    if (!ids) fail("Insufficient memory to validate IDs");
    int offset=0, k=0;
    for (int s=0; s<3; s++) {
        for (int i=0; i<a[s].n; i++) {
            if (a[s].p[i].id>INT_MAX-offset) fail("Particle ID overflow");
            a[s].p[i].id+=offset;
            ids[k++]=a[s].p[i].id;
        }
        offset+=a[s].n;
    }
    qsort(ids,(size_t)total,sizeof(int),compare_ids);
    for (int i=1; i<total; i++) if (ids[i]==ids[i-1]) fail("Duplicate merged particle ID");
    free(ids);
    /* Exclusive creation: never overwrite an input or previous intermediate. */
    FILE *out=fopen(argv[3],"wx");
    if (!out) { perror(argv[3]); return 2; }
    for (int t=0; t<4; t++) fprintf(out,"%d%c",a[0].counts[t]+a[1].counts[t]+a[2].counts[t],t==3?'\n':' ');
    for (int t=0; t<4; t++) for (int s=0; s<3; s++) for (int i=0; i<a[s].n; i++) {
        Particle *p=&a[s].p[i];
        if (p->type!=t) continue;
        fprintf(out,s==0 ? "%11d" : "%d",p->id);
        for (int j=0; j<8; j++) fprintf(out," %13.7g",(double)p->v[j]);
        fprintf(out," %d\n",p->type);
    }
    int bad=ferror(out);
    if (fclose(out)!=0 || bad) fail("Failed writing merged ASCII");
    for (int s=0; s<3; s++) free(a[s].p);
    return 0;
}
