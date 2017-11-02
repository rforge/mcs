#ifndef MCS_DETAIL_INTERRUPT_HH
#define MCS_DETAIL_INTERRUPT_HH



#ifndef MCS_CHECK_INTERRUPT
#define MCS_INTERRUPT_FLAG()   ((void) 0)
#define MCS_INTERRUPT_BREAK()  ((void) 0)
#else
#define MCS_INTERRUPT_FLAG()   MCS_CHECK_INTERRUPT()
#define MCS_INTERRUPT_BREAK()  if (MCS_INTERRUPT_FLAG())  break
#endif



#endif
